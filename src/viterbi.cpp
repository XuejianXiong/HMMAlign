#include "hmm_model.hpp"
#include <algorithm>
#include <limits>
#include <pybind11/pybind11.h>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace hmmalign;

struct AlignResult {
    float score;
    std::string path;
};

AlignResult viterbi_align_buffered(const std::string& read, const std::string& ref, 
                                  const ModelParams& params, AlignmentBuffer& buf) {
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const float neg_inf = -std::numeric_limits<float>::infinity();
    const size_t bandwidth = static_cast<size_t>(params.bandwidth);

    // 1. Allocate/Resize buffer (Capacity management)
    buf.reserve(total);

    // 2. Glocal Initialization
    // Reset only the first row to ensure clean start
    std::fill(buf.M.begin(), buf.M.begin() + width, neg_inf);
    std::fill(buf.I.begin(), buf.I.begin() + width, neg_inf);
    std::fill(buf.D.begin(), buf.D.begin() + width, neg_inf);
    
    for (size_t j = 0; j <= ref_len; ++j) {
        buf.M[j] = 0.0f; // Row 0: M can start anywhere for Glocal
    }

    // 3. Main DP Fill with Banding and Firewalls
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t j_start = (i > bandwidth ? i - bandwidth : 0);
        const size_t j_end = std::min(ref_len, i + bandwidth);
        const char char_i = read[i - 1];
        const size_t row_idx = i * width;
        const size_t prev_row_idx = (i - 1) * width;

        // FIREWALL: Clear the cell immediately before the band start
        // This prevents "stale" scores from previous buffer uses from leaking in.
        if (j_start > 0) {
            buf.M[row_idx + j_start - 1] = neg_inf;
            buf.I[row_idx + j_start - 1] = neg_inf;
            buf.D[row_idx + j_start - 1] = neg_inf;
        }

        for (size_t j = j_start; j <= j_end; ++j) {
            const size_t cur = row_idx + j;

            // Match State (Diagonal)
            if (j > 0) {
                const size_t diag = prev_row_idx + (j - 1);
                const float emit = (char_i == ref[j - 1]) ? emission.match_score : emission.mismatch_score;
                
                const float sM = buf.M[diag] + params.M_to_M;
                const float sI = buf.I[diag] + params.I_to_M;
                const float sD = buf.D[diag] + params.D_to_M;

                float best = sM;
                int8_t p = 0;
                if (sI > best) { best = sI; p = 1; }
                if (sD > best) { best = sD; p = 2; }

                if (best > neg_inf) {
                    buf.M[cur] = best + emit;
                    buf.ptrM[cur] = p;
                }
            }

            // Insertion (Up)
            const size_t up = prev_row_idx + j;
            const float sM_i = buf.M[up] + params.M_to_I;
            const float sI_i = buf.I[up] + params.I_to_I;
            
            buf.I[cur] = (sM_i >= sI_i) ? sM_i : sI_i;
            buf.ptrI[cur] = (sM_i >= sI_i) ? 0 : 1;

            // Deletion (Left)
            if (j > 0) {
                const size_t left = row_idx + (j - 1);
                const float sM_d = buf.M[left] + params.M_to_D;
                const float sD_d = buf.D[left] + params.D_to_D;
                
                buf.D[cur] = (sM_d >= sD_d) ? sM_d : sD_d;
                buf.ptrD[cur] = (sM_d >= sD_d) ? 0 : 2;
            }
        }
        
        // FIREWALL: Clear the cell immediately after the band end
        if (j_end < ref_len) {
            buf.M[row_idx + j_end + 1] = neg_inf;
            buf.I[row_idx + j_end + 1] = neg_inf;
            buf.D[row_idx + j_end + 1] = neg_inf;
        }
    }

    // 4. Termination: Search only the calculated band on the last row
    float final_score = neg_inf;
    size_t best_j = 0;
    int8_t state = -1;

    const size_t i_final = read_len;
    const size_t js = (i_final > bandwidth ? i_final - bandwidth : 0);
    const size_t je = std::min(ref_len, i_final + bandwidth);

    for (size_t j = js; j <= je; ++j) {
        size_t idx = i_final * width + j;
        if (buf.M[idx] > final_score) { final_score = buf.M[idx]; best_j = j; state = 0; }
        if (buf.I[idx] > final_score) { final_score = buf.I[idx]; best_j = j; state = 1; }
        if (buf.D[idx] > final_score) { final_score = buf.D[idx]; best_j = j; state = 2; }
    }
    
    if (final_score <= neg_inf || state == -1) return {neg_inf, ""};

    // 5. Traceback
    size_t ci = read_len, cj = best_j;
    std::string raw_path = "";
    while (ci > 0) {
        size_t cur = ci * width + cj;
        if (state == 0) { // Match
            raw_path += 'M';
            state = buf.ptrM[cur];
            ci--; cj--;
        } else if (state == 1) { // Insertion
            raw_path += 'I';
            state = buf.ptrI[cur];
            ci--;
        } else if (state == 2) { // Deletion
            raw_path += 'D';
            state = buf.ptrD[cur];
            cj--;
        } else {
            break;
        }
    }
    std::reverse(raw_path.begin(), raw_path.end());

    // 6. CIGAR Compression
    if (raw_path.empty()) return {final_score, ""};
    std::string cigar = "";
    char last_op = raw_path[0];
    int count = 0;
    for (char op : raw_path) {
        if (op == last_op) count++;
        else {
            cigar += std::to_string(count) + last_op;
            last_op = op;
            count = 1;
        }
    }
    cigar += std::to_string(count) + last_op;

    return {final_score, cigar};
}

// Simple wrapper for non-buffered calls
AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    AlignmentBuffer temp_buf;
    return viterbi_align_buffered(read, ref, params, temp_buf);
}

PYBIND11_MODULE(_core, m) {
    py::class_<ModelParams>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("M_to_M", &ModelParams::M_to_M)
        .def_readwrite("M_to_I", &ModelParams::M_to_I)
        .def_readwrite("M_to_D", &ModelParams::M_to_D)
        .def_readwrite("I_to_M", &ModelParams::I_to_M)
        .def_readwrite("I_to_I", &ModelParams::I_to_I)
        .def_readwrite("D_to_M", &ModelParams::D_to_M)
        .def_readwrite("D_to_D", &ModelParams::D_to_D)
        .def_readwrite("bandwidth", &ModelParams::bandwidth);

    py::class_<AlignmentBuffer>(m, "AlignmentBuffer").def(py::init<>());
    py::class_<AlignResult>(m, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("path", &AlignResult::path);

    m.def("viterbi_align_buffered", &viterbi_align_buffered);
    m.def("viterbi_align", &viterbi_align);
}
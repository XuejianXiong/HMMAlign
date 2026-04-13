#include <iostream>
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

/**
 * @brief Performs a Viterbi Alignment with true Soft-Clipping and Banding.
 * Includes a Local-Entry reset to handle leading/trailing adapter noise.
 */
AlignResult viterbi_align_buffered(const std::string& read, const std::string& ref, 
                                  const ModelParams& params, AlignmentBuffer& buf) {
    
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const float neg_inf = -std::numeric_limits<float>::infinity();
    const size_t bandwidth = static_cast<size_t>(params.bandwidth);

    buf.reserve(total);

    // 1. Initialization: Row 0 (Free start in reference)
    for (size_t j = 0; j <= ref_len; ++j) {
        buf.M[j] = 0.0f; 
        buf.I[j] = neg_inf;
        buf.D[j] = neg_inf;
    }

    // 2. Main DP Loop
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t row_idx = i * width;
        const size_t prev_row_idx = (i - 1) * width;
        const char char_i = read[i - 1];

        // Column 0: Free start in read for soft-clipping
        buf.M[row_idx] = 0.0f; 
        buf.I[row_idx] = neg_inf; 
        buf.D[row_idx] = neg_inf;

        const size_t j_start = (i > bandwidth ? i - bandwidth : 1);
        const size_t j_end = std::min(ref_len, i + bandwidth);

        // LEFT FIREWALL
        if (j_start > 0) {
            buf.M[row_idx + j_start - 1] = neg_inf;
            buf.I[row_idx + j_start - 1] = neg_inf;
            buf.D[row_idx + j_start - 1] = neg_inf;
        }

        for (size_t j = j_start; j <= j_end; ++j) {
            const size_t cur = row_idx + j;
            const char char_j = ref[j - 1];

            // --- Match State with Local Entry Reset ---
            const size_t diag = prev_row_idx + (j - 1);
            float emit = (char_i == 'N' || char_j == 'N') ? 0.0f : 
                         (char_i == char_j ? emission.match_score : emission.mismatch_score);
            
            float sM = buf.M[diag] + params.M_to_M;
            float sI = buf.I[diag] + params.I_to_M;
            float sD = buf.D[diag] + params.D_to_M;

            // Determine if starting fresh (0.0) is better than continuing path
            float best_m = std::max({sM, sI, sD, 0.0f}); 
            
            int8_t p_m = 0; 
            if (best_m <= 0.0f) {
                best_m = 0.0f; p_m = 3; // Source 3: New Alignment Boundary
            } else if (best_m == sI) {
                p_m = 1;
            } else if (best_m == sD) {
                p_m = 2;
            }

            buf.M[cur] = best_m + emit;
            buf.ptrM[cur] = p_m;

            // --- Insertion State ---
            const size_t up = prev_row_idx + j;
            float open_i = buf.M[up] + params.M_to_I;
            float extend_i = buf.I[up] + params.I_to_I;
            if (open_i >= extend_i) { buf.I[cur] = open_i; buf.ptrI[cur] = 0; }
            else { buf.I[cur] = extend_i; buf.ptrI[cur] = 1; }

            // --- Deletion State ---
            const size_t left = row_idx + (j - 1);
            float open_d = buf.M[left] + params.M_to_D;
            float extend_d = buf.D[left] + params.D_to_D;
            if (open_d >= extend_d) { buf.D[cur] = open_d; buf.ptrD[cur] = 0; }
            else { buf.D[cur] = extend_d; buf.ptrD[cur] = 2; }
        }
        
        // RIGHT FIREWALL
        if (j_end < ref_len) {
            buf.M[row_idx + j_end + 1] = neg_inf;
            buf.I[row_idx + j_end + 1] = neg_inf;
            buf.D[row_idx + j_end + 1] = neg_inf;
        }
    }

    // 3. Termination: Search for global max (Free-exit)
    float final_score = neg_inf;
    size_t best_i = 0, best_j = 0;
    int8_t state = -1;

    for (size_t i = 0; i <= read_len; ++i) {
        const size_t row_idx = i * width;
        const size_t js = (i > bandwidth ? i - bandwidth : 0);
        const size_t je = std::min(ref_len, i + bandwidth);
        for (size_t j = js; j <= je; ++j) {
            size_t idx = row_idx + j;
            if (buf.M[idx] > final_score) { 
                final_score = buf.M[idx]; 
                best_i = i; 
                best_j = j; 
                state = 0; 
            }
        }
    }
    
    if (state == -1) return {neg_inf, ""};

    // 4. Traceback
    size_t ci = best_i, cj = best_j;
    std::string raw_path = "";

    // Suffix Soft-Clipping
    if (best_i < read_len) {
        for (size_t k = 0; k < (read_len - best_i); ++k) raw_path += 'S';
    }

    while (ci > 0 && cj > 0) {
        size_t cur = ci * width + cj;
        if (state == 0) {
            int8_t next_state = buf.ptrM[cur];
            raw_path += 'M';
            ci--; cj--; // Consume the match
            if (next_state == 3) break; // Alignment start boundary reached
            state = next_state;
        }
        else if (state == 1) { raw_path += 'I'; state = buf.ptrI[cur]; ci--; }
        else if (state == 2) { raw_path += 'D'; state = buf.ptrD[cur]; cj--; }
    }

    // Prefix Soft-Clipping (remaining read bases)
    if (ci > 0) {
        for (size_t k = 0; k < ci; ++k) raw_path += 'S';
    }

    std::reverse(raw_path.begin(), raw_path.end());

    // 5. CIGAR Compression
    if (raw_path.empty()) return {final_score, ""};
    std::string cigar = "";
    char last_op = raw_path[0];
    int count = 0;
    for (char op : raw_path) {
        if (op == last_op) count++;
        else {
            cigar += std::to_string(count) + last_op;
            last_op = op; count = 1;
        }
    }
    cigar += std::to_string(count) + last_op;

    return {final_score, cigar};
}

/** Python Wrapper **/
AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    AlignmentBuffer temp_buf;
    return viterbi_align_buffered(read, ref, params, temp_buf);
}

/** Pybind11 Module **/
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
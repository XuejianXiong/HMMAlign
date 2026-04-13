#include "hmm_model.hpp"
#include <algorithm>
#include <limits>
#include <pybind11/pybind11.h>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace hmmalign;

static inline size_t index_2d(size_t i, size_t j, size_t width) {
    return i * width + j;
}

struct AlignResult {
    double score;
    std::string path;
};

// Core logic using the reused buffer
AlignResult viterbi_align_buffered(const std::string& read, const std::string& ref, 
                                  const ModelParams& params, AlignmentBuffer& buf) {
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const size_t bandwidth = params.bandwidth < 0 ? 0 : static_cast<size_t>(params.bandwidth);

    // Prepare buffer (only reallocates if necessary)
    buf.reserve(total);

    // --- 1. Glocal Initialization ---
    for (size_t j = 0; j <= ref_len; ++j) {
        buf.M[index_2d(0, j, width)] = 0.0;
    }

    // --- 2. Main DP Fill ---
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t j_start = (i > bandwidth ? i - bandwidth : 0);
        const size_t j_end = std::min(ref_len, i + bandwidth);

        for (size_t j = j_start; j <= j_end; ++j) {
            size_t cur = index_2d(i, j, width);

            if (j > 0) {
                size_t diag = index_2d(i - 1, j - 1, width);
                double emit = (read[i - 1] == ref[j - 1]) ? emission.match_score : emission.mismatch_score;
                double sM = buf.M[diag] + params.M_to_M;
                double sI = buf.I[diag] + params.I_to_M;
                double sD = buf.D[diag] + params.D_to_M;
                double best = std::max({sM, sI, sD});
                if (best > neg_inf) {
                    buf.M[cur] = best + emit;
                    if (sM >= sI && sM >= sD) buf.ptrM[cur] = 0;
                    else if (sI >= sD)         buf.ptrM[cur] = 1;
                    else                      buf.ptrM[cur] = 2;
                }
            }

            size_t up = index_2d(i - 1, j, width);
            double sM_i = buf.M[up] + params.M_to_I;
            double sI_i = buf.I[up] + params.I_to_I;
            double best_i = std::max(sM_i, sI_i);
            if (best_i > neg_inf) {
                buf.I[cur] = best_i;
                buf.ptrI[cur] = (sM_i >= sI_i) ? 0 : 1;
            }

            if (j > 0) {
                size_t left = index_2d(i, j - 1, width);
                double sM_d = buf.M[left] + params.M_to_D;
                double sD_d = buf.D[left] + params.D_to_D;
                double best_d = std::max(sM_d, sD_d);
                if (best_d > neg_inf) {
                    buf.D[cur] = best_d;
                    buf.ptrD[cur] = (sM_d >= sD_d) ? 0 : 2;
                }
            }
        }
    }

    // --- 3. Glocal Termination ---
    double final_score = neg_inf;
    size_t best_j = 0;
    int8_t state = -1;

    for (size_t j = 0; j <= ref_len; ++j) {
        size_t idx = index_2d(read_len, j, width);
        if (buf.M[idx] > final_score) { final_score = buf.M[idx]; best_j = j; state = 0; }
        if (buf.I[idx] > final_score) { final_score = buf.I[idx]; best_j = j; state = 1; }
        if (buf.D[idx] > final_score) { final_score = buf.D[idx]; best_j = j; state = 2; }
    }
    
    if (final_score == neg_inf) return {neg_inf, ""};

    // --- 4. Traceback ---
    size_t ci = read_len, cj = best_j;
    std::string raw_path = "";
    while (ci > 0) {
        size_t cur = index_2d(ci, cj, width);
        if (state == 0) { raw_path += 'M'; state = buf.ptrM[cur]; ci--; cj--; }
        else if (state == 1) { raw_path += 'I'; state = buf.ptrI[cur]; ci--; }
        else if (state == 2) { raw_path += 'D'; state = buf.ptrD[cur]; cj--; }
        else break;
        if (ci == 0) break;
    }
    std::reverse(raw_path.begin(), raw_path.end());

    // --- 5. CIGAR Compression ---
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

// Wrapper for the original function call (creates a temporary buffer)
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

    py::class_<AlignmentBuffer>(m, "AlignmentBuffer")
        .def(py::init<>());

    py::class_<AlignResult>(m, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("path", &AlignResult::path);

    m.def("viterbi_align", &viterbi_align, "Original version (allocates each time)");
    m.def("viterbi_align_buffered", &viterbi_align_buffered, "Fast version with buffer reuse");
}
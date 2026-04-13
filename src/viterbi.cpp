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

AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    const EmissionParams emission; // match_score = 1.0, mismatch = -1.5
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const size_t bandwidth = params.bandwidth < 0 ? 0 : static_cast<size_t>(params.bandwidth);

    std::vector<double> M(total, neg_inf), I(total, neg_inf), D(total, neg_inf);
    std::vector<int8_t> ptrM(total, -1), ptrI(total, -1), ptrD(total, -1);

    // --- 1. Glocal Initialization ---
    // Read can start anywhere in the reference (i=0) with 0.0 penalty.
    for (size_t j = 0; j <= ref_len; ++j) {
        M[index_2d(0, j, width)] = 0.0;
    }

    // --- 2. Main DP Fill ---
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t j_start = (i > bandwidth ? i - bandwidth : 0);
        const size_t j_end = std::min(ref_len, i + bandwidth);

        for (size_t j = j_start; j <= j_end; ++j) {
            size_t cur = index_2d(i, j, width);

            // Match/Mismatch State
            if (j > 0) {
                size_t diag = index_2d(i - 1, j - 1, width);
                double emit = (read[i - 1] == ref[j - 1]) ? emission.match_score : emission.mismatch_score;
                
                double sM = M[diag] + params.M_to_M;
                double sI = I[diag] + params.I_to_M;
                double sD = D[diag] + params.D_to_M;

                double best = std::max({sM, sI, sD});
                if (best > neg_inf) {
                    M[cur] = best + emit;
                    if (sM >= sI && sM >= sD) ptrM[cur] = 0;
                    else if (sI >= sD)         ptrM[cur] = 1;
                    else                      ptrM[cur] = 2;
                }
            }

            // Insertion (Gap in Reference)
            size_t up = index_2d(i - 1, j, width);
            double sM_i = M[up] + params.M_to_I;
            double sI_i = I[up] + params.I_to_I;
            double best_i = std::max(sM_i, sI_i);
            if (best_i > neg_inf) {
                I[cur] = best_i;
                ptrI[cur] = (sM_i >= sI_i) ? 0 : 1;
            }

            // Deletion (Gap in Read)
            if (j > 0) {
                size_t left = index_2d(i, j - 1, width);
                double sM_d = M[left] + params.M_to_D;
                double sD_d = D[left] + params.D_to_D;
                double best_d = std::max(sM_d, sD_d);
                if (best_d > neg_inf) {
                    D[cur] = best_d;
                    ptrD[cur] = (sM_d >= sD_d) ? 0 : 2;
                }
            }
        }
    }

    // --- 3. Glocal Termination: Find best score on last row ---
    double final_score = neg_inf;
    size_t best_j = 0;
    int8_t state = -1;

    for (size_t j = 0; j <= ref_len; ++j) {
        size_t idx = index_2d(read_len, j, width);
        if (M[idx] > final_score) { final_score = M[idx]; best_j = j; state = 0; }
        if (I[idx] > final_score) { final_score = I[idx]; best_j = j; state = 1; }
        if (D[idx] > final_score) { final_score = D[idx]; best_j = j; state = 2; }
    }
    
    if (final_score == neg_inf) return {neg_inf, ""};

    // --- 4. Traceback ---
    size_t ci = read_len, cj = best_j;
    std::string raw_path = "";
    
    while (ci > 0) {
        size_t cur = index_2d(ci, cj, width);
        int8_t current_state = state;

        if (state == 0) { // Match
            raw_path += 'M';
            state = ptrM[cur];
            ci--; cj--;
        } else if (state == 1) { // Insertion
            raw_path += 'I';
            state = ptrI[cur];
            ci--;
        } else if (state == 2) { // Deletion
            raw_path += 'D';
            state = ptrD[cur];
            cj--;
        } else {
            break;
        }
        if (ci == 0) break;
    }
    std::reverse(raw_path.begin(), raw_path.end());

    // --- 5. CIGAR Compression ---
    if (raw_path.empty()) return {final_score, ""};
    
    std::string cigar = "";
    char last_op = raw_path[0];
    int count = 0;
    for (char op : raw_path) {
        if (op == last_op) {
            count++;
        } else {
            cigar += std::to_string(count) + last_op;
            last_op = op;
            count = 1;
        }
    }
    cigar += std::to_string(count) + last_op;

    return {final_score, cigar};
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

    py::class_<AlignResult>(m, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("path", &AlignResult::path);

    m.def("viterbi_align", &viterbi_align, "Compute score and CIGAR path",
          py::arg("read"), py::arg("ref"), py::arg("params"));
}
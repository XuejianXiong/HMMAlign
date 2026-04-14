/**
 * @file viterbi.cpp
 * @brief High-performance implementation of the Viterbi algorithm for 
 * genomic alignment using an Affine Gap HMM.
 * * Architecture:
 * 1. Kernel: Banded DP recursion for score calculation.
 * 2. Traceback: Path reconstruction from pointer matrices.
 * 3. Compression: RAW path to CIGAR conversion.
 * 4. Bindings: Pybind11 interface.
 */

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

// --- Modular Helper Functions ---

/**
 * @brief Pre-calculates a score matrix for O(1) emission lookups.
 */
void init_score_matrix(float matrix[256][256], const EmissionParams& ep) {
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) {
            char a = std::toupper(static_cast<char>(i));
            char b = std::toupper(static_cast<char>(j));
            if (a == 'N' || b == 'N') matrix[i][j] = 0.0f;
            else if (a == b && (a == 'A' || a == 'C' || a == 'G' || a == 'T')) matrix[i][j] = ep.match_score;
            else matrix[i][j] = ep.mismatch_score;
        }
    }
}

/**
 * @brief Converts a raw path string (e.g., "MMMDD") into compressed CIGAR (3M2D).
 */
std::string compress_cigar(const std::string& raw_path) {
    if (raw_path.empty()) return "";
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
    return cigar;
}

/**
 * @brief Reconstructs the optimal alignment path.
 */
std::string perform_traceback(size_t best_i, size_t best_j, int8_t state, 
                             size_t read_len, size_t width, const AlignmentBuffer& buf) {
    size_t ci = best_i, cj = best_j;
    std::string raw_path = "";

    // End-of-read soft-clipping
    if (best_i < read_len) raw_path.append(read_len - best_i, 'S');

    while (ci > 0 && cj > 0) {
        size_t cur = ci * width + cj;
        if (state == 0) { // Match State
            int8_t next_state = buf.ptrM[cur];
            raw_path += 'M';
            ci--; cj--;
            if (next_state == 3) break; // Reached local start
            state = next_state;
        }
        else if (state == 1) { raw_path += 'I'; state = buf.ptrI[cur]; ci--; }
        else if (state == 2) { raw_path += 'D'; state = buf.ptrD[cur]; cj--; }
    }

    // Start-of-read soft-clipping
    if (ci > 0) raw_path.append(ci, 'S');
    std::reverse(raw_path.begin(), raw_path.end());
    return raw_path;
}

// --- Main Engine ---

AlignResult viterbi_align_buffered(const std::string& read, const std::string& ref, 
                                  const ModelParams& params, AlignmentBuffer& buf) {
    
    const EmissionParams emission_params;
    static float score_matrix[256][256];
    static bool matrix_init = false;
    if (!matrix_init) {
        init_score_matrix(score_matrix, emission_params);
        matrix_init = true;
    }

    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const float neg_inf = -std::numeric_limits<float>::infinity();
    const size_t bandwidth = static_cast<size_t>(params.bandwidth);

    buf.reserve((read_len + 1) * width);

    // 1. Initialization
    for (size_t j = 0; j <= ref_len; ++j) {
        buf.M[j] = 0.0f; buf.I[j] = neg_inf; buf.D[j] = neg_inf;
    }

    // 2. DP Recursion
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t row_idx = i * width;
        const size_t prev_row_idx = (i - 1) * width;
        const uint8_t char_i = static_cast<uint8_t>(read[i - 1]);

        buf.M[row_idx] = 0.0f; // Soft-clip init
        const size_t j_start = (i > bandwidth ? i - bandwidth : 1);
        const size_t j_end = std::min(ref_len, i + bandwidth);

        for (size_t j = j_start; j <= j_end; ++j) {
            const size_t cur = row_idx + j;
            const size_t diag = prev_row_idx + (j - 1);
            const size_t up = prev_row_idx + j;
            const size_t left = row_idx + (j - 1);
            
            // Match State Logic
            float emit = score_matrix[char_i][static_cast<uint8_t>(ref[j - 1])];
            float sM = buf.M[diag] + params.M_to_M;
            float sI = buf.I[diag] + params.I_to_M;
            float sD = buf.D[diag] + params.D_to_M;

            float best_m = std::max({sM, sI, sD, 0.0f}); 
            int8_t p_m = (best_m <= 0.0f) ? 3 : (best_m == sI ? 1 : (best_m == sD ? 2 : 0));

            buf.M[cur] = (best_m <= 0.0f ? 0.0f : best_m) + emit;
            buf.ptrM[cur] = p_m;

            // Affine Insertion
            float open_i = buf.M[up] + params.M_to_I;
            float ext_i = buf.I[up] + params.I_to_I;
            if (open_i >= ext_i) { buf.I[cur] = open_i; buf.ptrI[cur] = 0; }
            else { buf.I[cur] = ext_i; buf.ptrI[cur] = 1; }

            // Affine Deletion
            float open_d = buf.M[left] + params.M_to_D;
            float ext_d = buf.D[left] + params.D_to_D;
            if (open_d >= ext_d) { buf.D[cur] = open_d; buf.ptrD[cur] = 0; }
            else { buf.D[cur] = ext_d; buf.ptrD[cur] = 2; }
        }
    }

    // 3. Termination & Path Reconstruction
    float final_score = neg_inf;
    size_t best_i = 0, best_j = 0;
    for (size_t i = 0; i <= read_len; ++i) {
        const size_t r_idx = i * width;
        const size_t js = (i > bandwidth ? i - bandwidth : 0);
        const size_t je = std::min(ref_len, i + bandwidth);
        for (size_t j = js; j <= je; ++j) {
            if (buf.M[r_idx + j] > final_score) { 
                final_score = buf.M[r_idx + j]; best_i = i; best_j = j; 
            }
        }
    }
    
    if (final_score == neg_inf) return {neg_inf, ""};
    
    std::string raw = perform_traceback(best_i, best_j, 0, read_len, width, buf);
    return {final_score, compress_cigar(raw)};
}

/** Python Wrapper **/
AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    AlignmentBuffer temp_buf;
    return viterbi_align_buffered(read, ref, params, temp_buf);
}

/** Pybind11 Module Definition **/
PYBIND11_MODULE(_core, m) {
    m.doc() = "HMMAlign Core: Optimized Viterbi engine.";
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

    py::class_<AlignmentBuffer>(m, "AlignmentBuffer").def(py::init<>()).def("reserve", &AlignmentBuffer::reserve);
    py::class_<AlignResult>(m, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("path", &AlignResult::path);

    m.def("viterbi_align_buffered", &viterbi_align_buffered);
    m.def("viterbi_align", &viterbi_align);
}
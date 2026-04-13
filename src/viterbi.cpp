#include <iostream>
#include "hmm_model.hpp"
#include <algorithm>
#include <limits>
#include <pybind11/pybind11.h>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace hmmalign;

/**
 * @brief Result structure for Python binding
 */
struct AlignResult {
    float score;
    std::string path;
};

/**
 * @brief Performs a Glocal Viterbi Alignment with Affine Gap Penalties and Banding.
 * * This implementation uses a pre-allocated buffer to minimize heap churn and
 * incorporates "Firewalls" (boundary sanitization) to allow safe buffer reuse.
 */
AlignResult viterbi_align_buffered(const std::string& read, const std::string& ref, 
                                  const ModelParams& params, AlignmentBuffer& buf) {
    
    std::cout << "DEBUG C++: Bandwidth received = " << params.bandwidth << std::endl;
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const float neg_inf = -std::numeric_limits<float>::infinity();
    // const size_t bandwidth = static_cast<size_t>(params.bandwidth);
    const size_t bandwidth = 200; // FORCED HARD-CODE

    // 1. Memory Management: Ensure buffer is large enough for current read/ref pair
    buf.reserve(total);

    // 2. Initialization: Set up Row 0 for Glocal alignment
    // M[0][j] = 0 allows the alignment to start at any position in the reference without penalty.
    for (size_t j = 0; j <= ref_len; ++j) {
        buf.M[j] = 0.0f; 
        buf.I[j] = neg_inf;
        buf.D[j] = neg_inf;
    }

    // 3. Main Dynamic Programming Loop
    for (size_t i = 1; i <= read_len; ++i) {
        const size_t row_idx = i * width;
        const size_t prev_row_idx = (i - 1) * width;
        const char char_i = read[i - 1];

        // Column 0 Sanitization: Prevents alignment from starting mid-read at ref index 0
        buf.M[row_idx] = neg_inf;
        buf.I[row_idx] = neg_inf; 
        buf.D[row_idx] = neg_inf;

        // Determine Banding Bounds: O(N * Bandwidth) complexity
        const size_t j_start = (i > bandwidth ? i - bandwidth : 1);
        const size_t j_end = std::min(ref_len, i + bandwidth);

        // LEFT FIREWALL: Explicitly kill stale data from previous buffer usage
        if (j_start > 0) {
            buf.M[row_idx + j_start - 1] = neg_inf;
            buf.I[row_idx + j_start - 1] = neg_inf;
            buf.D[row_idx + j_start - 1] = neg_inf;
        }

        for (size_t j = j_start; j <= j_end; ++j) {
            const size_t cur = row_idx + j;
            const char char_j = ref[j - 1];

            // --- Match State Transition ---
            const size_t diag = prev_row_idx + (j - 1);
            
            // N-Neutrality: Ambiguous bases provide 0 reward/penalty to preserve path sensitivity
            float emit;
            if (char_i == 'N' || char_j == 'N') {
                emit = 0.0f;
            } else {
                emit = (char_i == char_j) ? emission.match_score : emission.mismatch_score;
            }
            
            const float sM = buf.M[diag] + params.M_to_M;
            const float sI = buf.I[diag] + params.I_to_M;
            const float sD = buf.D[diag] + params.D_to_M;

            float best_m = sM; int8_t p_m = 0;
            if (sI > best_m) { best_m = sI; p_m = 1; }
            if (sD > best_m) { best_m = sD; p_m = 2; }

            if (best_m > neg_inf) {
                buf.M[cur] = best_m + emit;
                buf.ptrM[cur] = p_m;
            }

            // --- Insertion State (Gap in Reference) ---
            const size_t up = prev_row_idx + j;
            const float open_i = buf.M[up] + params.M_to_I;
            const float extend_i = buf.I[up] + params.I_to_I;
            
            if (open_i >= extend_i) {
                buf.I[cur] = open_i; buf.ptrI[cur] = 0; // From Match (Gap Open)
            } else {
                buf.I[cur] = extend_i; buf.ptrI[cur] = 1; // From Insertion (Gap Extend)
            }

            // --- Deletion State (Gap in Read) ---
            const size_t left = row_idx + (j - 1);
            const float open_d = buf.M[left] + params.M_to_D;
            const float extend_d = buf.D[left] + params.D_to_D;
            
            if (open_d >= extend_d) {
                buf.D[cur] = open_d; buf.ptrD[cur] = 0; // From Match (Gap Open)
            } else {
                buf.D[cur] = extend_d; buf.ptrD[cur] = 2; // From Deletion (Gap Extend)
            }
        }
        
        // RIGHT FIREWALL: Ensure the band doesn't bleed into stale memory
        if (j_end < ref_len) {
            buf.M[row_idx + j_end + 1] = neg_inf;
            buf.I[row_idx + j_end + 1] = neg_inf;
            buf.D[row_idx + j_end + 1] = neg_inf;
        }
    }

    // 4. Termination: Find the best score in the final row's valid band
    float final_score = neg_inf;
    size_t best_j = 0;
    int8_t state = -1;
    const size_t i_f = read_len;
    const size_t js = (i_f > bandwidth ? i_f - bandwidth : 0);
    const size_t je = std::min(ref_len, i_f + bandwidth);

    for (size_t j = js; j <= je; ++j) {
        size_t idx = i_f * width + j;
        if (buf.M[idx] > final_score) { final_score = buf.M[idx]; best_j = j; state = 0; }
        if (buf.I[idx] > final_score) { final_score = buf.I[idx]; best_j = j; state = 1; }
        if (buf.D[idx] > final_score) { final_score = buf.D[idx]; best_j = j; state = 2; }
    }
    
    if (final_score <= neg_inf || state == -1) return {neg_inf, ""};

    // 5. Traceback: Reconstruct path from end to start
    size_t ci = read_len, cj = best_j;
    std::string raw_path = "";
    while (ci > 0 || (ci == 0 && state == 2)) {
        size_t cur = ci * width + cj;
        if (state == 0) { raw_path += 'M'; state = buf.ptrM[cur]; ci--; cj--; }
        else if (state == 1) { raw_path += 'I'; state = buf.ptrI[cur]; ci--; }
        else if (state == 2) { raw_path += 'D'; state = buf.ptrD[cur]; cj--; }
        else break;
        if (ci == 0 && state != 2) break; // Terminate if row 0 reached (unless deleting)
    }
    std::reverse(raw_path.begin(), raw_path.end());

    // 6. CIGAR Compression: Convert raw path (MMMDMM) to CIGAR (3M1D2M)
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

/**
 * @brief Python-friendly wrapper for non-buffered alignment calls
 */
AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    AlignmentBuffer temp_buf;
    return viterbi_align_buffered(read, ref, params, temp_buf);
}

/**
 * @brief Pybind11 Module Definition
 */
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

    m.def("viterbi_align_buffered", &viterbi_align_buffered, "Buffered Viterbi alignment");
    m.def("viterbi_align", &viterbi_align, "Standard Viterbi alignment");
}
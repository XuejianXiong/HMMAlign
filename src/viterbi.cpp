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

// Result structure to return to Python
struct AlignResult {
    double score;
    std::string path;
};

AlignResult viterbi_align(std::string read, std::string ref, ModelParams params) {
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const double neg_inf = -std::numeric_limits<double>::infinity();

    std::vector<double> M(total, neg_inf), I(total, neg_inf), D(total, neg_inf);
    // 0: Match, 1: Insert, 2: Delete
    std::vector<int8_t> ptrM(total, -1), ptrI(total, -1), ptrD(total, -1);

    M[0] = 0.0;

    for (size_t i = 0; i <= read_len; ++i) {
        for (size_t j = 0; j <= ref_len; ++j) {
            if (i == 0 && j == 0) continue;
            size_t cur = index_2d(i, j, width);

            // --- Match State ---
            if (i > 0 && j > 0) {
                size_t diag = index_2d(i - 1, j - 1, width);
                double emit = (read[i - 1] == ref[j - 1]) ? emission.match_score : emission.mismatch_score;
                
                double sM = M[diag] + params.M_to_M;
                double sI = I[diag] + params.I_to_M;
                double sD = D[diag] + params.D_to_M;

                if (sM >= sI && sM >= sD) { M[cur] = sM + emit; ptrM[cur] = 0; }
                else if (sI >= sD)         { M[cur] = sI + emit; ptrM[cur] = 1; }
                else                      { M[cur] = sD + emit; ptrM[cur] = 2; }
            }

            // --- Insertion State (Gap in Ref) ---
            if (i > 0) {
                size_t up = index_2d(i - 1, j, width);
                double sM = M[up] + params.M_to_I;
                double sI = I[up] + params.I_to_I;

                if (sM >= sI) { I[cur] = sM; ptrI[cur] = 0; }
                else          { I[cur] = sI; ptrI[cur] = 1; }
            }

            // --- Deletion State (Gap in Read) ---
            if (j > 0) {
                size_t left = index_2d(i, j - 1, width);
                double sM = M[left] + params.M_to_D;
                double sD = D[left] + params.D_to_D;

                if (sM >= sD) { D[cur] = sM; ptrD[cur] = 0; }
                else          { D[cur] = sD; ptrD[cur] = 2; }
            }
        }
    }

    // --- Traceback ---
    size_t ci = read_len, cj = ref_len;
    std::string path = "";
    
    // Pick starting state
    size_t last = index_2d(ci, cj, width);
    int8_t state; 
    double final_score;
    if (M[last] >= I[last] && M[last] >= D[last]) { state = 0; final_score = M[last]; }
    else if (I[last] >= D[last])                 { state = 1; final_score = I[last]; }
    else                                         { state = 2; final_score = D[last]; }

    while (ci > 0 || cj > 0) {
        size_t cur = index_2d(ci, cj, width);
        if (state == 0) { path += 'M'; state = ptrM[cur]; ci--; cj--; }
        else if (state == 1) { path += 'I'; state = ptrI[cur]; ci--; }
        else if (state == 2) { path += 'D'; state = ptrD[cur]; cj--; }
    }
    std::reverse(path.begin(), path.end());

    return {final_score, path};
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
        .def_readwrite("D_to_D", &ModelParams::D_to_D);

    py::class_<AlignResult>(m, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("path", &AlignResult::path);

    m.def("viterbi_align", &viterbi_align, "Compute score and path",
          py::arg("read"), py::arg("ref"), py::arg("params"));
}
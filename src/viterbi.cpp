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

static double max3(double a, double b, double c) {
    return std::max(a, std::max(b, c));
}

double viterbi_align(std::string read, std::string ref, ModelParams params) {
    const EmissionParams emission;
    const size_t read_len = read.size();
    const size_t ref_len = ref.size();
    const size_t width = ref_len + 1;
    const size_t total = (read_len + 1) * width;
    const double neg_inf = -std::numeric_limits<double>::infinity();

    std::vector<double> M(total, neg_inf);
    std::vector<double> I(total, neg_inf);
    std::vector<double> D(total, neg_inf);

    // Initialization: only the origin may be reached with score 0.
    M[0] = 0.0;

    for (size_t i = 1; i <= read_len; ++i) {
        for (size_t j = 1; j <= ref_len; ++j) {
            const size_t cur = index_2d(i, j, width);
            const size_t diag = index_2d(i - 1, j - 1, width);
            const size_t up = index_2d(i - 1, j, width);
            const size_t left = index_2d(i, j - 1, width);

            double emit = (read[i - 1] == ref[j - 1]) ? emission.match_score : emission.mismatch_score;
            double m_score = max3(
                M[diag] + params.M_to_M,
                I[diag] + params.I_to_M,
                D[diag] + params.D_to_M);
            M[cur] = m_score + emit;

            double i_score = std::max(
                M[up] + params.M_to_I,
                I[up] + params.I_to_I);
            I[cur] = i_score;

            double d_score = std::max(
                M[left] + params.M_to_D,
                D[left] + params.D_to_D);
            D[cur] = d_score;
        }
    }

    const size_t final_index = index_2d(read_len, ref_len, width);
    return std::max({M[final_index], I[final_index], D[final_index]});
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "HMMAlign C++ Core Engine";
    py::class_<ModelParams>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("M_to_M", &ModelParams::M_to_M)
        .def_readwrite("M_to_I", &ModelParams::M_to_I)
        .def_readwrite("M_to_D", &ModelParams::M_to_D)
        .def_readwrite("I_to_M", &ModelParams::I_to_M)
        .def_readwrite("I_to_I", &ModelParams::I_to_I)
        .def_readwrite("D_to_M", &ModelParams::D_to_M)
        .def_readwrite("D_to_D", &ModelParams::D_to_D);

    m.def("viterbi_align", &viterbi_align,
          "Compute the Viterbi alignment score for a read against a reference",
          py::arg("read"), py::arg("ref"), py::arg("params"));
}

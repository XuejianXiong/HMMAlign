#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

// Placeholder function to test the bridge
double align_test(std::string read, std::string ref) {
    // We will implement the real Viterbi logic here next
    return (double)(read.length() + ref.length());
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "HMMAlign C++ Core Engine";
    m.def("align_test", &align_test, "A placeholder function to test the C++ bridge",
          py::arg("read"), py::arg("ref"));
}
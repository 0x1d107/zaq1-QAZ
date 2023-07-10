#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include "inv_q_filter.cpp"

// PYBIND11_MODULE(inv_q_filtr, m) {
//     m.doc() = "Inverse Q filtration";
//     m.def("phase_q_correction", &phase_q_correction);
// }
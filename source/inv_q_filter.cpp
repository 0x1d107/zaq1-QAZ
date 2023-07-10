// #include <pybind11/pybind11.h>
// #include <pybind11/eigen.h>
// #include <pybind11/stl.h>
// #include <pybind11/complex.h>
// #include <pybind11/functional.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "inv_q_filter.hpp"
#include "stft.hpp"

std::vector<double> phase_q_correction(const std::vector<double> &signal, double fs, double Q, double target_freq, int nperseg, int noverlap) {
    int signal_size = signal.size();

    std::vector<double> freq_samples = get_freq_samples(nperseg, fs);

    std::vector<double> time_samples = get_time_samples(signal_size, nperseg, noverlap, fs);

    int time_samples_count = time_samples.size();

    auto original_tf_spectra = stft(signal, nperseg, noverlap);

    Eigen::MatrixXcd corrected_tf_spectra(nperseg, time_samples_count);
    
    const std::complex<double> I (0.0, 1.0);

    double gamma = 1. / Q / M_PI;
    // ЧТобы сигнал не менялся

    // for (int i = 0; i < nperseg / 2 + 1; ++i) {
    //     for (int j = 0; j < time_samples_count; ++j) {
    //         std::complex<double> correction_term ((std::pow(((double) freq_samples[i] / target_freq), gamma) - 1) * 2. * M_PI * freq_samples[i] * time_samples[j], 0.);
    //         corrected_tf_spectra(i, j) = original_tf_spectra(i, j); 
    //     }
    // }

    for (int j = 0; j < time_samples_count; ++j) {
        corrected_tf_spectra(0, j) = I;
    }

    for (int i = 1; i < nperseg / 2 + 1; ++i) {
        for (int j = 0; j < time_samples_count; ++j) {
           std::complex<double> correction_term = std::exp(I * (std::pow((freq_samples[i] /  target_freq), -gamma) - 1.0) * 2.0 * M_PI * freq_samples[i] * time_samples[j]);
           corrected_tf_spectra(i, j) = correction_term * original_tf_spectra(i, j); 
        }
    }

    return istft(corrected_tf_spectra, nperseg, noverlap);
}


// PYBIND11_MODULE(inv_q_filtr, m) {
//     m.doc() = "Inverse Q filtration";
//     m.def("phase_q_correction", &phase_q_correction);
// }
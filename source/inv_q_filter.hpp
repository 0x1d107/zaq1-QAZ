#pragma once

#include <vector>

std::vector<double> phase_q_correction(const std::vector<double>& signal, double fs, double Q, double target_freq, int nperseg, int noverlap);

// std::vector<double> amp_q_correction(const std::vector<double>& signal, double fs, double Q, double target_freq, int nperseg, int noverlap);
// std::vector<double> amp_q_correction(const std::vector<double>& signal, double fs, double Q, double sigma2, double target_freq, int nperseg, int noverlap);
// std::vector<double> amp_q_correction(const std::vector<double>& signal, double fs, double Q, double gain_limit, double target_freq, int nperseg, int noverlap);
// std::vector<double> amp_q_correction(const std::vector<double>& signal, double fs, double Q, double sigma2, double gain_limit, double target_freq, int nperseg, int noverlap);

// std::vector<double> full_q_correction  (const std::vector<double>& signal, double fs, double Q, double target_freq, int nperseg, int noverlap);
// std::vector<double> full_q_correction  (const std::vector<double>& signal, double fs, double Q, double sigma2, double target_freq, int nperseg, int noverlap);
// std::vector<double> full_q_correction  (const std::vector<double>& signal, double fs, double Q, double gain_limit, double target_freq, int nperseg, int noverlap);
// std::vector<double> full_q_correction  (const std::vector<double>& signal, double fs, double Q, double sigma2, double gain_limit, double target_freq, int nperseg, int noverlap);
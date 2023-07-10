#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/FFT>
#include <complex>

typedef double (*wfunc_t)(int n,int N);
typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> CM;

extern wfunc_t window;

CM stft(const std::vector<double> &signal,int nperseg = 256,int overlap = 128);
std::vector<double>  istft(const CM& spectrum,int nperseg = 256,int overlap = 128);
void getTimeFreq(double fs,int nperseg, int overlap,int len, std::vector<double> &time, std::vector<double> &frq);

std::vector<double> get_freq_samples(int nperseg, int fs);
std::vector<double>  get_time_samples(int signal_length, int nperseg, int noverlap, double fs);

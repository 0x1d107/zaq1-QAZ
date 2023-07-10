// #include <pybind11/pybind11.h>
// #include <pybind11/eigen.h>
// #include <pybind11/stl.h>
// #include <pybind11/complex.h>
// #include <pybind11/functional.h>

#include <eigen3/unsupported/Eigen/FFT>
#include <complex>

#include "stft.hpp"

double hann_window(int n,int N){
	if(n<0||n>N)
		return 0;
	double s = sin(M_PI*n/(N-1));
	return s*s;
}
double fullspec_window(int n,int N){
	return 1;
}
wfunc_t window = hann_window;
CM stft(const std::vector<double> &signal,int nperseg,int overlap){
	
	//
	Eigen:Eigen::FFT<double>fft;
	int window_step = nperseg - overlap;
	int T = (signal.size()/window_step)+((signal.size()%window_step)>0)+1;
	CM stftm(T,nperseg);
	for(int i=0;i<T;i++){
		std::vector<double> f(nperseg);
		for(int j=0;j<nperseg;j++){
			f[j] = window(j,nperseg) * signal[j+i*window_step];
		}
		std::vector<std::complex<double>> frq(nperseg,0);
		fft.fwd(frq,f);
		for(int j=0;j<nperseg;j++){
			stftm(i,j) = frq[j];
		}


	}
	stftm.transposeInPlace();
	return stftm;
}

std::vector<double>  istft(const CM& spectrum,int nperseg,int noverlap){
	CM transposed_spectrum(spectrum);
	transposed_spectrum.transposeInPlace();
	Eigen::FFT<double> fft;
	int window_step = nperseg - noverlap;
	int T = transposed_spectrum.rows();
	int length = (T-1)*window_step+nperseg;
	std::vector<double> signal(length);
	for(int i=0;i<T;i++){
		std::vector<std::complex<double>> frq(nperseg,0);
		std::vector<double> f(nperseg);
		for(int j = 0; j < nperseg;j++)
			frq[j] = transposed_spectrum(i,j);
		fft.inv(f,frq);
		for(int j = 0; j < nperseg;j++)
			if(window(j,nperseg)>0)signal[j+window_step*i] = f[j]/window(j,nperseg);
	
	}

	return signal;	
}

int get_freq_count(int nperseg) {
	if (nperseg % 2 == 0) {
		return nperseg / 2 + 1;
	}
	else {
		return (nperseg - 1) / 2 + 1;
	}
}


std::vector<double> get_freq_samples(int nperseg, int fs) {
	std::vector<double> freq (nperseg, 0.);

	for (int i = 0; i < nperseg / 2 + 1; ++i) {
		freq[i] = ((double) i / nperseg * fs);
	}

	for (int i = nperseg - 1, j = 1; i > nperseg / 2; --i, ++j) {
		freq[i] = -((double) j / nperseg * fs);
	}
	
	return freq;
}

int  get_time_count(int signal_length, int nperseg, int noverlap) {
	return std::ceil((double) signal_length / (nperseg - noverlap)) + 1;
}

std::vector<double>  get_time_samples(int signal_length, int nperseg, int noverlap, double fs) {
	int times_count = get_time_count(signal_length, nperseg, noverlap);

	std::vector<double>  time_samples (times_count);
	
	for (int i = 0; i < time_samples.size(); ++i) {
		time_samples[i] = i * noverlap / fs;
	}

	return time_samples;
}

// PYBIND11_MODULE(inv_q_filtr, m) {
//   m.doc() = "Short time fourier transform";

//   m.def("stft", &stft);

//   m.def("istft", &istft);

//   m.def("get_freq_samples", &get_freq_samples);

//   m.def("get_time_samples", &get_time_samples);
// }
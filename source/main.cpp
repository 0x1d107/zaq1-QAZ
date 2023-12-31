#include <iostream>
#include <fstream>

#include "inv_q_filter.hpp"


int main() {

    std::ifstream myfile;
	myfile.open("signal.txt");
	std::vector<double> signal;
	if (myfile.is_open()) { 
		double val;
		while(myfile >> val)  signal.push_back(val);
	} else {
		std::cout << "NONONONO\n";
	}

    auto corrected = phase_q_correction(signal, 1000, 200, 20, 256, 128);
	std::ofstream myfile2;
	myfile2.open ("result.txt");
	for (int i = 0; i < corrected.size(); ++i) {
		myfile2 << corrected[i] << std::endl;
	}
	myfile2.close();
}

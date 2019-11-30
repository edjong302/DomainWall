#include "time_evolution.hpp"
#include <iostream>
using namespace std;

std::vector<float> time_evolution::slope(std::vector<float> i_vector){
	std::vector<float> result = i_vector;
	return result;
}

/*std::vector<float> time_evolution::slope_rk_1_psi(std::vector<float> i_vector, float i_dx){
	std::vector<float> result;
	int size = i_vector.size();
	std::cout << size << std::endl;
	for (int i = 0; i < size; i++){
		result[i] = (i_vector[(i - 1 + size)%size] - 2*i_vector[i] + i_vector[(i + 1)%size])/(i_dx*i_dx) - potential.dPotentialdfield(i_vector[i]);
		//result[i] = (1.)/(i_dx*i_dx);
	}
	std::cout << size << std::endl;
		result [i] = 0.0;
	
	return result;
}*/
#include "time_evolution.hpp"
#include <iostream>
#include <vector>
using namespace std;

std::vector<float> time_evolution::slope_phi_rk1(std::vector<float> i_vector){
	std::vector<float> result = i_vector;
	return result;
}

std::vector<float> time_evolution::slope_psi_rk1(std::vector<float> i_vector, float i_dx){
	int size = i_vector.size();
	std::vector<float> result(size);
	for (int i = 0; i < size - 1; i++){
		result[i] = (i_vector[(i - 1 + size)] - 2*i_vector[i] + i_vector[(i + 1)])/(i_dx*i_dx);
	}
	return result;
}
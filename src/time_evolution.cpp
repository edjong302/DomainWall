#include "time_evolution.hpp"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

time_evolution::time_evolution(){
	lambda = 2000.0;
	eta = 1.3;
	woof = 250.0;
}

float time_evolution::the_Potential(float field){
	return lambda/4.0*pow((pow(field,2.) - pow(eta,2.)),2) - woof/(2.*eta)*(field - eta);
}

float time_evolution::dPotentialdfield(float field){
	return lambda*field*(pow(field,2) - pow(eta,2)) - woof/(2*eta)*field;
}

float time_evolution::get_lambda(){
	return lambda;
}

float time_evolution::get_eta(){
	return eta;
}

std::vector<float> time_evolution::rk1_phi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	return l0_phi(i_vector, ii_vector, i_dx, i_alpha);
}

std::vector<float> time_evolution::rk1_psi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	return l0_psi(i_vector, ii_vector, i_dx, i_alpha);
}

std::vector<float> time_evolution::rk2_phi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	int size = i_vector.size();
	std::vector<float> result(size);
	for (int i = 0; i < size; i++){
			result[i] = (l0_phi(i_vector, ii_vector, i_dx, i_alpha)[i] + l1_phi(i_vector, ii_vector, i_dx, i_alpha)[i])/2;
		}
	return result;
}

std::vector<float> time_evolution::rk2_psi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	int size = i_vector.size();
	std::vector<float> result(size);
	for (int i = 0; i < size; i++){
			result[i] = (l0_psi(i_vector, ii_vector, i_dx, i_alpha)[i] + l1_psi(i_vector, ii_vector, i_dx, i_alpha)[i])/2;
		}
	return result;
}

std::vector<float> time_evolution::l0_phi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	std::vector<float> result = ii_vector;
	return result;
}

std::vector<float> time_evolution::l0_psi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	int size = i_vector.size();
	std::vector<float> result(size);
	for (int i = 0; i < size; i++){
		result[i] = (i_vector[(i - 1 + size)%size] - 2*i_vector[i] + i_vector[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(i_vector[i]);
	}
	// 	result[i] = (i_vector[(i - 1 + size)%size] - 2*i_vector[i] + i_vector[(i + 1)%size])/(i_dx*i_dx) - dPotentialdfield(i_vector[i]);
	return result;
}

std::vector<float> time_evolution::l1_phi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	int size = i_vector.size();
	std::vector<float> result(size);
	for (int i = 0; i < size - 1; i++){
		result[i] = ii_vector[i] + i_alpha*i_dx*l0_psi(i_vector, ii_vector, i_dx, i_alpha)[i];
		}
	return result;
}

std::vector<float> time_evolution::l1_psi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha){
	int size = i_vector.size();
	std::vector<float> result(size);
	std::vector<float> phi_temp(size);
	for (int i = 0; i < size - 1; i++){
		phi_temp[i] = i_vector[i] + i_alpha*i_dx*l0_phi(i_vector, ii_vector, i_dx, i_alpha)[i];
		result[i] = (phi_temp[(i - 1 + size)%size] - 2*phi_temp[i] + phi_temp[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(phi_temp[i]);
		}
	return result;
}
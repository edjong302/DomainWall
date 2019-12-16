#include "time_evolution.hpp"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

time_evolution::time_evolution(){
	
}

double time_evolution::the_Potential(double field){ // This is not actually important for the evolution.
	return lambda/4.0*pow((pow(field,2.) - pow(eta,2.)),2) - woof/(2.*eta)*(field - eta);
}

double time_evolution::dPotentialdfield(double field){ // This runs the evolution of the field.
	return lambda*pow(field,3) + woof*pow(field, 2) + eta*field;
}

double time_evolution::get_lambda(){
	return lambda;
}

double time_evolution::get_eta(){
	return eta;
}

void time_evolution::set_lambda(double i_lambda){
	lambda = i_lambda;
}

void time_evolution::set_eta(double i_eta){
	eta = i_eta;
}

void time_evolution::set_woof(double i_woof){
	woof = i_woof;
}

void time_evolution::set_kreiss(double i_kreiss){
	kreiss = i_kreiss;
}

std::vector<double> time_evolution::rk1_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	return l0_phi(i_vector, ii_vector, i_dx, i_alpha);
}

std::vector<double> time_evolution::rk1_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	return l0_psi(i_vector, ii_vector, i_dx, i_alpha);
}

std::vector<double> time_evolution::rk2_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> l0_v = l0_phi(i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l1_v = l1_phi(l0_v, i_vector, ii_vector, i_dx, i_alpha);
	for (int i = 0; i < size; i++){
		//result[i] = (l0_v[i] + l1_v[i])/2;
		result[i] = l1_v[i];
	}
	return result;
}

std::vector<double> time_evolution::rk2_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> l0_v = l0_psi(i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l1_v = l1_psi(l0_v, i_vector, ii_vector, i_dx, i_alpha);
	for (int i = 0; i < size; i++){
		//result[i] = (l0_v[i] + l1_v[i])/2;
		result[i] = l1_v[i];
	}
	return result;
}

std::vector<double> time_evolution::rk4_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> l0_v = l0_phi(i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l1_v = l1_phi(l0_v, i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l2_v = l2_phi(l1_v, i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l3_v = l3_phi(l2_v, i_vector, ii_vector, i_dx, i_alpha);
	for (int i = 0; i < size; i++){
		//result[i] = (l0_v[i] + l1_v[i])/2;
		result[i] = (l0_v[i] + 2.*l1_v[i] + 2.*l2_v[i] + l3_v[i])/6.;
		result[i] += kreiss/(64.*i_dx)*(i_vector[(i + 3)%size] - 6.*i_vector[(i + 2)%size] + 15.*i_vector[(i + 1)%size] -20.*i_vector[i] + 15.*i_vector[(i - 1 + size)%size] - 6.*i_vector[(i - 2 + size)%size] + i_vector[(i - 3 + size)%size]);
	}
	return result;
}

std::vector<double> time_evolution::rk4_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> l0_v = l0_psi(i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l1_v = l1_psi(l0_v, i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l2_v = l2_psi(l1_v, i_vector, ii_vector, i_dx, i_alpha);
	std::vector<double> l3_v = l3_psi(l2_v, i_vector, ii_vector, i_dx, i_alpha);
	for (int i = 0; i < size; i++){
		//result[i] = (l0_v[i] + l1_v[i])/2;
		result[i] = (l0_v[i] + 2.*l1_v[i] + 2.*l2_v[i] + l3_v[i])/6.;
		result[i] += kreiss/(64.*i_dx)*(ii_vector[(i + 3)%size] - 6.*ii_vector[(i + 2)%size] + 15.*ii_vector[(i + 1)%size] -20.*ii_vector[i] + 15.*ii_vector[(i - 1 + size)%size] - 6.*ii_vector[(i - 2 + size)%size] + ii_vector[(i - 3 + size)%size]);
	}
	return result;
}

std::vector<double> time_evolution::l0_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	std::vector<double> result = ii_vector;
	return result;
}

std::vector<double> time_evolution::l0_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	for (int i = 0; i < size; i++){
		result[i] = (i_vector[(i + size - 1)%size] - 2*i_vector[i] + i_vector[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(i_vector[i]);
	}
	return result;
}

std::vector<double> time_evolution::l1_phi(std::vector<double> i_l0, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	for (int i = 0; i < size; i++){
		result[i] = ii_vector[i] + i_alpha*i_dx/2.*i_l0[i];
		}
	return result;
}

std::vector<double> time_evolution::l1_psi(std::vector<double> i_l0, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> phi_temp(size);
	for (int i = 0; i < size; i++){
		phi_temp[i] = i_vector[i] + i_alpha*i_dx/2.*i_l0[i];
	}
	for (int i = 0; i < size; i++){
		result[i] = (phi_temp[(i + size - 1)%size] - 2*phi_temp[i] + phi_temp[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(phi_temp[i]);
	}
	return result;
}

std::vector<double> time_evolution::l2_phi(std::vector<double> i_l1, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	for (int i = 0; i < size; i++){
		result[i] = ii_vector[i] + i_alpha*i_dx/2.*i_l1[i];
		}
	return result;
}

std::vector<double> time_evolution::l2_psi(std::vector<double> i_l1, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> phi_temp(size);
	for (int i = 0; i < size; i++){
		phi_temp[i] = i_vector[i] + i_alpha*i_dx/2.*i_l1[i];
	}
	for (int i = 0; i < size; i++){
		result[i] = (phi_temp[(i + size - 1)%size] - 2*phi_temp[i] + phi_temp[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(phi_temp[i]);
	}
	return result;
}

std::vector<double> time_evolution::l3_phi(std::vector<double> i_l2, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	for (int i = 0; i < size; i++){
		result[i] = ii_vector[i] + i_alpha*i_dx*i_l2[i];
		}
	return result;
}

std::vector<double> time_evolution::l3_psi(std::vector<double> i_l2, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha){
	int size = i_vector.size();
	std::vector<double> result(size);
	std::vector<double> phi_temp(size);
	for (int i = 0; i < size; i++){
		phi_temp[i] = i_vector[i] + i_alpha*i_dx*i_l2[i];
	}
	for (int i = 0; i < size; i++){
		result[i] = (phi_temp[(i + size - 1)%size] - 2*phi_temp[i] + phi_temp[(i + 1)%size])/(i_dx*i_dx) -1.*dPotentialdfield(phi_temp[i]);
	}
	return result;
}
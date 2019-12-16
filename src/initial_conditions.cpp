#include "initial_conditions.hpp"
#include <math.h>

double initial_conditions::get_offset1(){
	return offset1;
}

double initial_conditions::get_offset2(){
	return offset2;
}

double initial_conditions::get_sigma(){
	return sigma;
}

void initial_conditions::set_offset1(double i_offset1){
	offset1 = i_offset1;
}

void initial_conditions::set_offset2(double i_offset2){
	offset2 = i_offset2;
}

void initial_conditions::set_sigma(double i_sigma){
	sigma = i_sigma;
}

double initial_conditions::phi_init(double i_offset1, double i_offset2, double i_sigma, double i_x){
	// Domain wall:
	double vacuum_field_value = 2.01981;
	return vacuum_field_value*.5*(-tanh((i_x - i_offset1)/i_sigma) + tanh((i_x - i_offset2)/i_sigma) + 2);
	//return -0.01*(tanh((i_x - i_offset1)/i_sigma) - tanh((i_x - i_offset2)/i_sigma) - 1);
	// Vacuum with small perturbation:
	//return 0.01 + 0.0001*exp(-pow(-(i_x - 100.),2.)/100.);
	// Stable vacuum:
	//return 0.01;
	// Finetuned for asymmetric potential:
	// double min1 = -0.009304;
	// double min2 = 0.01057;
	// return min2 + 0.5*(min1 - min2)*tanh((i_x - i_offset1)/i_sigma) - 0.5*(min1 - min2)*tanh((i_x - i_offset2)/i_sigma);
}

double initial_conditions::psi_init(double i_offset1, double i_offset2, double i_sigma, double i_x){
	return 0.;
}
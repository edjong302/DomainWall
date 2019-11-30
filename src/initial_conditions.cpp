#include "initial_conditions.hpp"
#include <math.h>

float initial_conditions::get_offset1(){
	return offset1;
}

float initial_conditions::get_offset2(){
	return offset2;
}

float initial_conditions::get_sigma(){
	return sigma;
}

float initial_conditions::phi_init(float i_offset1, float i_offset2, float i_sigma, float i_x){
	return -0.01*(tanh((i_x - i_offset1)/i_sigma) - tanh((i_x - i_offset2)/i_sigma) - 1);
	//return 0.0101;
}

float initial_conditions::psi_init(float i_offset1, float i_offset2, float i_sigma, float i_x){
	return 0.;
}
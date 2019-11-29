#include "potential.hpp"

Potential::Potential(){
	lambda = 2000.0;
	eta = 0.01;
	epsilon = 0.0001;
}

float Potential::the_Potential(float field){
	return lambda/4.0*(field*field - eta*eta)*(field*field - eta*eta) - epsilon/(2*eta)*(field - eta);
}

float Potential::dPotentialdfield(float field){
	return lambda*field*(field*field - eta*eta) - epsilon/(2*eta);
}

float Potential::get_lambda(){
	return lambda;
}

float Potential::get_eta(){
	return eta;
}
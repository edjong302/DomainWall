#ifndef POTENTIAL
#define POTENTIAL

class Potential{
	private:
		float lambda, eta, epsilon;
	public:
		Potential();
		float get_lambda();
		float get_eta();
		float the_Potential(float field);
		float dPotentialdfield(float field);
};

#endif
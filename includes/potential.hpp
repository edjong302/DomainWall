#ifndef POTENTIAL
#define POTENTIAL

class Potential{
	private:
		float lambda, eta, woof;
	public:
		Potential();
		float get_lambda();
		float get_eta();
		float the_Potential(float field);
		float dPotentialdfield(float field);
};

#endif
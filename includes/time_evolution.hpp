#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H
#include <vector>
using namespace std;


class time_evolution{
	private:
		float lambda, eta, woof;
	public:
		time_evolution();
		float get_lambda();
		float get_eta();
		float the_Potential(float field);
		float dPotentialdfield(float field);
		std::vector<float> rk1_phi(std::vector<float> i_vector);
		std::vector<float> rk1_psi(std::vector<float> i_vector, float i_dx);
		std::vector<float> l0_phi(std::vector<float> i_vector);
		std::vector<float> l0_psi(std::vector<float> i_vector, float i_dx);
		std::vector<float> l1_phi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha);
		std::vector<float> l1_psi(std::vector<float> i_vector, std::vector<float> ii_vector, float i_dx, float i_alpha);
};

#endif
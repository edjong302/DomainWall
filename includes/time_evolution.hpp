#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H
#include <vector>
using namespace std;


class time_evolution{
	private:
	public:
		std::vector<float> slope_phi_rk1(std::vector<float> i_vector);
		std::vector<float> slope_psi_rk1(std::vector<float> i_vector, float i_dx);
};

#endif
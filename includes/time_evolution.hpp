#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H
#include <vector>
using namespace std;


class time_evolution{
	private:
		double lambda, eta, woof;
		int kreiss;
	public:
		time_evolution();
		double get_lambda();
		double get_eta();
		double the_Potential(double field);
		double dPotentialdfield(double field);
		std::vector<double> rk1_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> rk1_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> rk2_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> rk2_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> rk4_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> rk4_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l0_phi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l0_psi(std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l1_phi(std::vector<double> i_l0, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l1_psi(std::vector<double> i_l0, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l2_phi(std::vector<double> i_l1, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l2_psi(std::vector<double> i_l1, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l3_phi(std::vector<double> i_l2, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
		std::vector<double> l3_psi(std::vector<double> i_l2, std::vector<double> i_vector, std::vector<double> ii_vector, double i_dx, double i_alpha);
};

#endif
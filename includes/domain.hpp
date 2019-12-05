#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <string>
#include "initial_conditions.hpp"
#include "time_evolution.hpp"

class DomainWall{
    private:
        int L, N;
        double alpha, t_max;
        std::vector<double> new_state_phi;
        std::vector<double> current_state_phi;
        std::vector<double> new_state_psi;
        std::vector<double> current_state_psi;
        initial_conditions init_cond;
        time_evolution time_evol;
    public:
        DomainWall();
        int get_L();
        std::vector<double> get_new_state_phi();
        std::vector<double> get_current_state_phi();
        void set_initial_conditions();
        void take_step(double i_steps);
        void evolve_to_the_end();
        void write_phi(std::string filename);
        void write_grid(std::string filename);
        void write_time(std::string filename);
        void delete_files(std::string filename1, std::string filename2, std::string filename3);
};

#endif

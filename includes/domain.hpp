#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <string>
#include "initial_conditions.hpp"
#include "time_evolution.hpp"

class DomainWall{
    private:
        int L, N;
        float alpha, t_max;
        std::vector<float> new_state_phi;
        std::vector<float> current_state_phi;
        std::vector<float> new_state_psi;
        std::vector<float> current_state_psi;
        initial_conditions init_cond;
        time_evolution time_evol;
    public:
        DomainWall();
        int get_L();
        std::vector<float> get_new_state_phi();
        std::vector<float> get_current_state_phi();
        void set_initial_conditions();
        void take_step(float i_steps);
        void evolve_to_the_end();
        void write_phi(std::string filename);
        void write_grid(std::string filename);
        void write_time(std::string filename);
        void delete_files(std::string filename1, std::string filename2, std::string filename3);
};

#endif

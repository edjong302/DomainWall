#include "domain.hpp"
#include "initial_conditions.hpp"
#include "time_evolution.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
using namespace std;

DomainWall::DomainWall(){
    L = 200;
    N = 500;
    alpha = 0.1;
    t_max = 100;
    std::vector<float> ivec(N);
    std::vector<float> iivec(N);
    std::vector<float> iiivec(N);
    std::vector<float> ivvec(N);
    current_state_phi = ivec;
    new_state_phi = iivec;
    current_state_psi = iiivec;
    new_state_psi = ivvec;
};

int DomainWall::get_L(){
    return L;
};

std::vector<float> DomainWall::get_current_state_phi(){
    return current_state_phi;
}

std::vector<float> DomainWall::get_new_state_phi(){
    return new_state_phi;
}

void DomainWall::set_initial_conditions(){
    float dx = float(L)/N;
	float offset1 = 0.3*L;
	float offset2 = 0.7*L;
	float sigma = 0.02*L;
    for (int i = 0; i < N; i++){
        current_state_phi[i] = init_cond.phi_init(offset1, offset2, sigma, i*dx);
        new_state_phi[i] = init_cond.phi_init(offset1, offset2, sigma, i*dx);
        current_state_psi[i] = 0.0;
        new_state_psi[i] = 0.0;
    }
    write_phi("output.txt");
}

void DomainWall::evolve_to_the_end(){
    float max_steps = (t_max*N)/(alpha*L);
    std::cout << "max steps"<< " " << max_steps << std::endl;
    float steps = 0;
    while (steps < max_steps){
        take_step(steps);
        write_phi("output.txt");
        steps += 1.;
    }
    std::string Fliss = "woof";
    std::cout << Fliss << std::endl;
    
}

void DomainWall::take_step(float i_steps){
    current_state_phi = new_state_phi;
    current_state_psi = new_state_psi;
    float dx = float(L)/(N);
    std::vector<float> phi_slope = time_evol.rk2_phi(current_state_phi, current_state_psi, dx, alpha);
    std::vector<float> psi_slope = time_evol.rk2_psi(current_state_phi, current_state_psi, dx, alpha);
    for (int i = 0; i < N; i++){
        // if (i_steps == 4.){
        //     std::cout << psi_slope[i] << std::endl;
        // }
        new_state_phi[i] = current_state_phi[i] + alpha*dx*phi_slope[i];
        //new_state_psi[i] = current_state_psi[i] + alpha*(float(N) + 1.0)/L*(current_state_phi[(i + N - 1)%N] - 2*current_state_phi[i] + current_state_phi[(i + 1)%N]) - alpha*float(L)/(N + 1)*lambda*current_state_phi[i%N]*(current_state_phi[i]*current_state_phi[i] - eta*eta);
        new_state_psi[i] = current_state_psi[i] + alpha*dx*psi_slope[i];
    }
    //std::cout << i_steps << std::endl;
}

void DomainWall::delete_files(std::string filename1, std::string filename2, std::string filename3){
    std::ofstream outfile1;
    std::ofstream outfile2;
    std::ofstream outfile3;
    outfile1.open(filename1);
    outfile1.close();
    outfile2.open(filename2);
    outfile2.close();
    outfile3.open(filename3);
    outfile3.close();
}

void DomainWall::write_phi(std::string filename){
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    for (int i = 0; i < N; i++){
        outfile << current_state_phi[i] << " ";
    }
    outfile << "\n";
}

void DomainWall::write_grid(std::string filename){
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    for (int i = 0; i < N; i++){
        outfile << i*float(L)/(N + 1) << std::endl;
    }
    outfile.close();
}

void DomainWall::write_time(std::string filename){
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    float dx = float(L)/(N + 1);
    int N_max = t_max/(alpha*dx);
    float dt = alpha*dx;
    for (int i = 0; i < N_max + 1; i++){
        outfile << i*dt << std::endl;
    }
    outfile.close();
}

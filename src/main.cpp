#include "domain.hpp"
#include <iostream>
using namespace std;
#include<fstream>

int main(){
    DomainWall solution;
    solution.delete_files("grid.txt", "output.txt", "time.txt");
    solution.write_grid("grid.txt");
    std::cout << 1 << std::endl;
    solution.write_time("time.txt");
    std::cout << 1 << std::endl;
    solution.set_initial_conditions();
    std::cout << 1 << std::endl;
    solution.evolve_to_the_end();
    // for (int i = 0; i < 50; i++){
    //     std::cout << solution.get_current_state_phi()[i] << std::endl;
    // }
    // solution.take_step();
    // for (int i = 0; i < 50; i++){
    //     std::cout << solution.get_new_state_phi()[i] << std::endl;
    // }
    // solution.evolve_to_the_end();
    // for (int i = 0; i<50; i++){
    //     std::cout << solution.get_new_state_phi()[i] << std::endl;
    // }
    return 0;
}
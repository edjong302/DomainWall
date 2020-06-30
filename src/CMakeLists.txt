add_library(domain domain.cpp)
add_library(initial_conditions initial_conditions.cpp)
add_library(potential potential.cpp)
add_library(time_evolution time_evolution.cpp)

add_executable(main main.cpp)
target_link_libraries(main domain
                           initial_conditions
                           potential
                           time_evolution
                     )
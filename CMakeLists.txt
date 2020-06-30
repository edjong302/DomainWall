project(DomainWall)
cmake_minimum_required(VERSION 3.0)
include_directories(includes)

enable_testing()
add_subdirectory(src)
add_subdirectory(tests)
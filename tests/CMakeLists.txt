add_subdirectory(/usr/src/gtest ${PROJECT_BINARY_DIR}/gtests)

include(CTest)

add_executable(test_potential test_potential.cpp)
target_link_libraries(test_potential
                      gtest
                      gtest_main
                      potential)
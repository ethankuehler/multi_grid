cmake_minimum_required(VERSION 3.14)
project(c_SOR_3d C)

set(CMAKE_C_STANDARD 11)

string(APPEND CMAKE_C_FLAGS "-fopenmp")
add_subdirectory(multi_grid)
add_subdirectory(tSolver)
add_executable(c_SOR_3d main.c solve_block.c multi_grid/operators.c multi_grid/multi_grid.c)
target_link_libraries(c_SOR_3d m multi_grid tSolver)

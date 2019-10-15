cmake_minimum_required(VERSION 3.14)
project(c_SOR_3d C CXX)

set(CMAKE_C_STANDARD 11)

string(APPEND CMAKE_C_FLAGS_DEBUG "")
SET(CMAKE_C_FLAGS -pthread)

add_subdirectory(helpers)

add_executable(c_SOR_3d main.c solve_block.c multi_grid.c)
target_link_libraries(c_SOR_3d m solver multi_helpers)
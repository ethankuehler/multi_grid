add_library(multi_helpers STATIC operators.c N_len.c solve_thread.c)
add_library(solver solve_thread.c)

link_libraries(solver multi_helpers)
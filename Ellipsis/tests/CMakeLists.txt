cmake_minimum_required(VERSION 2.6)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -O3")
endif(CMAKE_COMPILER_IS_GNUCC)

add_executable(test_c_ellipsis.exe test_hanson.c test_kinetic_energy.c test_mt19937.c test_c_ellipsis.c )
target_link_libraries(test_c_ellipsis.exe ellipsis -lm)

add_test(guided_hmc test_c_ellipsis.exe)

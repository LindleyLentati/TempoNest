cmake_minimum_required(VERSION 2.6)

project(ellipsis_examples Fortran C)

#get the fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER}")

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
	set(CMAKE_COMPILER_IS_GNUFC true)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	set(CMAKE_COMPILER_IS_INTELFC true)
else(Fortran_COMPILER_NAME STREQUAL "gfortran")
	if(NOT FORCE_UNKNOWN_FORTRAN_COMPILER)
		message(FATAL_ERROR "Untested FORTRAN Compiler ${Fortran_COMPILER_NAME}")
	endif(NOT FORCE_UNKNOWN_FORTRAN_COMPILER)	
endif(Fortran_COMPILER_NAME STREQUAL "gfortran")


#set the C compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -O3")	
elseif(CMAKE_COMPILER_IS_INTELCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC ")	
endif(CMAKE_COMPILER_IS_GNUCC)

#set the fortran compiler flags
if(CMAKE_COMPILER_IS_GNUFC)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f95 -pedantic -Wall -Wextra  -Wfatal-errors -O3")
elseif(CMAKE_COMPILER_IS_INTELFC)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f90 -warn general -O3")
endif(CMAKE_COMPILER_IS_GNUFC)

add_executable(uncorr_gauss_c.exe uncorr_gauss_c.c)
target_link_libraries(uncorr_gauss_c.exe ellipsis -lm)

add_executable(uncorr_gauss_f.exe uncorr_gauss_f.f90)
target_link_libraries(uncorr_gauss_f.exe ellipsis)

#for the correlated gaussian example, we need LAPACK
find_package(LAPACK)
if(LAPACK_FOUND)
	add_executable(gauss_f.exe lin_alg.f90 gauss_f_mod.f90 gauss_f.f90)
	target_link_libraries(gauss_f.exe ellipsis ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

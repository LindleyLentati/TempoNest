cmake_minimum_required(VERSION 2.6)

project(ellipsis C )

if(ENABLE_TESTING)
	enable_testing()
endif(ENABLE_TESTING)

# guard against in-source builds (got this from Eigen)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

#see what compilers are we using
#get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)

#if(NOT CMAKE_COMPILER_IS_GNUCC)
#	if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
#		set(CMAKE_COMPILER_IS_INTELCC true)
#	else(${CMAKE_C_COMPILER} MATCHES "icc.*$")
#		if(NOT FORCE_UNKNOWN_C_COMPILER)
#			message(FATAL_ERROR "Untested C Compiler ${C_COMPILER_NAME}")
#		endif(NOT FORCE_UNKNOWN_C_COMPILER)
#	endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
#endif(NOT CMAKE_COMPILER_IS_GNUCC)


include_directories(ellipsis)
add_subdirectory(ellipsis)
add_subdirectory(examples)
if(ENABLE_TESTING)
	add_subdirectory(tests)
endif(ENABLE_TESTING)

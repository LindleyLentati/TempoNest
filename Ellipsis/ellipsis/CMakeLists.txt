cmake_minimum_required(VERSION 2.6)


#if(CMAKE_COMPILER_IS_GNUCC)
#	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC")
#elseif(CMAKE_COMPILER_IS_INTELCC)
#	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC ")
#endif(CMAKE_COMPILER_IS_GNUCC)

#set compiler flags
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c89 -pedantic -Wall -Wextra -Wfatal-errors -fPIC")
	include_directories(/usr/include/x86_64-linux-gnu/c++/4.7)#my ubuntu has a problem
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

add_library(ellipsis guided_hmc.c  hanson.c  mt19937.c)


# Guided Hamiltonian Sampler

[![Build Status](https://travis-ci.org/tbs1980/Ellipsis.svg?branch=master)](https://travis-ci.org/tbs1980/Ellipsis)
[![codecov.io](http://codecov.io/github/tbs1980/Ellipsis/coverage.svg?branch=master)](http://codecov.io/github/tbs1980/Ellipsis?branch=master)

## Authors

0. [Sreekumar Thaithara Balan](mailto:sbalan@star.ucl.ac.uk)
1. [Michael Hobson](mailto:mph@mrao.cam.ac.uk)
2. [Mark Ashdown](mailto:maja1.mrao.cam.ac.uk)

## Compilation

We need CMake to compile the files. Using git clone the repository

	$ git clone https://github.com/tbs1980/Ellipsis.git

Then move to the directory `Ellipsis` and follow the instructions below.

	$ cd Ellipsis
	$ mkdir build
	$ cd build
	$ cmake ../
	$ make

The library `libellipsis.a` will be in the directory `ellipsis` and the examples
can be found in `examples`.

	$ ls -la ellipsis
	total 40
	drwxrwxr-x 3 sbalan sbalan    84 Jul  1 12:06 .
	drwxrwxr-x 5 sbalan sbalan   115 Jul  1 12:06 ..
	drwxrwxr-x 3 sbalan sbalan    84 Jul  1 12:06 CMakeFiles
	-rw-rw-r-- 1 sbalan sbalan  1171 Jul  1 12:06 cmake_install.cmake
	-rw-rw-r-- 1 sbalan sbalan 26830 Jul  1 12:06 libellipsis.a
	-rw-rw-r-- 1 sbalan sbalan  7604 Jul  1 12:06 Makefile

	$ ls -la examples
	drwxrwxr-x 3 sbalan sbalan  4096 Jul  1 12:06 .
	drwxrwxr-x 5 sbalan sbalan   115 Jul  1 12:06 ..
	drwxrwxr-x 5 sbalan sbalan   145 Jul  1 12:06 CMakeFiles
	-rw-rw-r-- 1 sbalan sbalan  1171 Jul  1 12:06 cmake_install.cmake
	-rwxrwxr-x 1 sbalan sbalan 37888 Jul  1 12:06 gauss_f.exe
	-rw-rw-r-- 1 sbalan sbalan  8103 Jul  1 12:06 gauss_post.mod
	-rw-rw-r-- 1 sbalan sbalan  2097 Jul  1 12:06 linearalgebrautils.mod
	-rw-rw-r-- 1 sbalan sbalan  8723 Jul  1 12:06 Makefile
	-rwxrwxr-x 1 sbalan sbalan 25096 Jul  1 12:06 uncorr_gauss_c.exe
	-rwxrwxr-x 1 sbalan sbalan 27017 Jul  1 12:06 uncorr_gauss_f.exe

## LICENSE

The code is distributed under [MPL2](https://www.mozilla.org/MPL/2.0/)

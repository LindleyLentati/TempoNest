SET(CMAKE_Fortran_COMPILER "/usr/local/x86_64/gnu/gcc-4.8.2/bin/gfortran")
SET(CMAKE_Fortran_COMPILER_ARG1 "")
SET(CMAKE_Fortran_COMPILER_ID "GNU")
SET(CMAKE_Fortran_PLATFORM_ID "")

SET(CMAKE_AR "/usr/bin/ar")
SET(CMAKE_RANLIB "/usr/bin/ranlib")
SET(CMAKE_COMPILER_IS_GNUG77 1)
SET(CMAKE_Fortran_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

SET(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_Fortran_COMPILER_ID_RUN 1)
SET(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;f77;F77;f90;F90;for;For;FOR;f95;F95)
SET(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_Fortran_LINKER_PREFERENCE 20)
IF(UNIX)
  SET(CMAKE_Fortran_OUTPUT_EXTENSION .o)
ELSE(UNIX)
  SET(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
ENDIF(UNIX)

# Save compiler ABI information.
SET(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
SET(CMAKE_Fortran_COMPILER_ABI "")
SET(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

IF(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
ENDIF()

IF(CMAKE_Fortran_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
ENDIF(CMAKE_Fortran_COMPILER_ABI)

IF(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  SET(CMAKE_LIBRARY_ARCHITECTURE "")
ENDIF()

SET(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;quadmath;m;c")
SET(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/usr/local/cuda-6.0/lib64;/usr/local/x86_64/gnu/gcc-4.8.2/lib/gcc/x86_64-unknown-linux-gnu/4.8.2;/usr/local/x86_64/gnu/gcc-4.8.2/lib64;/lib64;/usr/lib64;/usr/local/intel/composer_xe_2011_sp1.7.256/compiler/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/ipp/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21;/usr/local/cuda-6.0/lib;/usr/local/x86_64/gnu/openmpi-1.8.3-dlopen/lib;/usr/local/x86_64/gnu/mpc-1.0/lib;/usr/local/x86_64/gnu/mpfr-3.1.1/lib;/usr/local/x86_64/gnu/gmp-5.0.5/lib;/usr/local/intel-15.3.0/composer_xe_2015.3.187/ipp/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/compiler/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/mkl/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/tbb/lib/intel64/gcc4.4;/usr/local/x86_64/gnu/cfitsio-3.290/lib;/usr/local/x86_64/gnu/gcc-4.8.2/lib")

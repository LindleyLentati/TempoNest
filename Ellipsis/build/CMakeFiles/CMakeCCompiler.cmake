SET(CMAKE_C_COMPILER "/usr/local/x86_64/gnu/gcc-4.8.2/bin/gcc")
SET(CMAKE_C_COMPILER_ARG1 "")
SET(CMAKE_C_COMPILER_ID "GNU")
SET(CMAKE_C_COMPILER_VERSION "4.8.2")
SET(CMAKE_C_PLATFORM_ID "Linux")

SET(CMAKE_AR "/usr/bin/ar")
SET(CMAKE_RANLIB "/usr/bin/ranlib")
SET(CMAKE_LINKER "/usr/bin/ld")
SET(CMAKE_COMPILER_IS_GNUCC 1)
SET(CMAKE_C_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_C_COMPILER_ENV_VAR "CC")

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_C_COMPILER_ID_RUN 1)
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
SET(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
SET(CMAKE_C_SIZEOF_DATA_PTR "8")
SET(CMAKE_C_COMPILER_ABI "ELF")
SET(CMAKE_C_LIBRARY_ARCHITECTURE "")

IF(CMAKE_C_SIZEOF_DATA_PTR)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
ENDIF(CMAKE_C_SIZEOF_DATA_PTR)

IF(CMAKE_C_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
ENDIF(CMAKE_C_COMPILER_ABI)

IF(CMAKE_C_LIBRARY_ARCHITECTURE)
  SET(CMAKE_LIBRARY_ARCHITECTURE "")
ENDIF()

SET(CMAKE_C_HAS_ISYSROOT "")


SET(CMAKE_C_IMPLICIT_LINK_LIBRARIES "pthread;c")
SET(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/usr/local/cuda-6.0/lib64;/usr/local/x86_64/gnu/gcc-4.8.2/lib/gcc/x86_64-unknown-linux-gnu/4.8.2;/usr/local/x86_64/gnu/gcc-4.8.2/lib64;/lib64;/usr/lib64;/usr/local/intel/composer_xe_2011_sp1.7.256/compiler/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/ipp/lib/intel64;/usr/local/intel/composer_xe_2011_sp1.7.256/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21;/usr/local/cuda-6.0/lib;/usr/local/x86_64/gnu/openmpi-1.8.3-dlopen/lib;/usr/local/x86_64/gnu/mpc-1.0/lib;/usr/local/x86_64/gnu/mpfr-3.1.1/lib;/usr/local/x86_64/gnu/gmp-5.0.5/lib;/usr/local/intel-15.3.0/composer_xe_2015.3.187/ipp/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/compiler/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/mkl/lib/intel64;/usr/local/intel-15.3.0/composer_xe_2015.3.187/tbb/lib/intel64/gcc4.4;/usr/local/x86_64/gnu/cfitsio-3.290/lib;/usr/local/x86_64/gnu/gcc-4.8.2/lib")

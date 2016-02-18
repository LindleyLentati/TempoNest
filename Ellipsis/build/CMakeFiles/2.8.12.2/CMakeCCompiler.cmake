set(CMAKE_C_COMPILER "/usr/local/Cluster-Apps/gcc/4.8.1/bin/gcc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "4.8.1")
set(CMAKE_C_PLATFORM_ID "Linux")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()




set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/usr/local/Cluster-Apps/cuda/6.0/lib64;/usr/local/Cluster-Apps/cuda/6.0/extras/CUPTI/lib64;/usr/local/Cluster-Apps/cuda/6.0/nvvm/lib64;/usr/local/Cluster-Apps/gcc/4.8.1/lib/gcc/x86_64-unknown-linux-gnu/4.8.1;/usr/local/Cluster-Apps/gcc/4.8.1/lib64;/lib64;/usr/lib64;/usr/local/Cluster-Apps/intel/cce/15.0.1.133/ipp/lib/intel64;/usr/local/Cluster-Apps/intel/cce/15.0.1.133/tbb/lib/intel64;/usr/local/Cluster-Apps/intel/fce/15.0.1.133/lib/intel64;/usr/local/Cluster-Apps/pgplot;/usr/local/Cluster-Apps/cuda/6.0/lib;/usr/local/Cluster-Apps/cuda/6.0/open64/lib;/usr/local/Cluster-Apps/fftw/intel/64/3.3.3/lib;/usr/local/Cluster-Apps/intel/mkl/10.3.10.319/composer_xe_2011_sp1.10.319/mkl/lib/intel64;/usr/local/Cluster-Apps/intel/impi/4.1.3.045/lib64;/usr/local/Cluster-Users/sjr20/cfitsio/3.03/icc/lib;/usr/local/Cluster-Apps/gcc/4.8.1/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")




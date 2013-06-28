#
# SWIN_LIB_CULA([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the NVIDIA cuda/cula GPU programming
# compilers and libraries.
#
# CULA_CFLAGS - autoconfig variable with flags required for compiling
# CULA_LIBS   - autoconfig variable with flags required for linking
# HAVE_CULA   - automake conditional
# HAVE_CULA   - pre-processor macro in config.h
#
# Some environment variables are required for the CUDA/CULA libraries to
# function, if they are not installed in standard locations
#
# CULA_LIB_PATH_64 or CULA_LIB_PATH_32 need to be pointing to the installation
# directory (64 bit will be tried first)
# CULA_INC_PATH needs to point to the headers
#
# This macro tries to link a test program in the following way
#
#    -L${CULA_LIB_PATH_64} -lcula_core -lcula_lapack -lcula_lapack_fortran -lcublas -lcudart -lcuda -I$(CULA_INC)
#
# TODO: Test whether the BLAS and LAPACK libraries really are required
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_CULA],
[
  AC_PROVIDE([SWIN_LIB_CULA])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])

  AC_MSG_CHECKING([for CUDA/CULA installation])

  CULA_CFLAGS=""
  CULA_LIBS=""

  CULA_LIB="-lcula_lapack -lcuda"

  if test x"$CULA_LIB_PATH_64" != x; then
    CULA_LIBS="-L$CULA_LIB_PATH_64"
  else
    if test x"$CULA_LIB_PATH_32" != x; then
      CULA_LIBS="-L$CULA_LIB_PATH_32"
    fi
  fi

  if test x"$CULA_INC_PATH" != x; then
    CULA_CFLAGS="-I$CULA_INC_PATH"
  fi

  CULA_LIBS="$CULA_LIBS $CULA_LIB"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $CULA_LIBS"
  CFLAGS="$ac_save_CFLAGS $CULA_CFLAGS"

  AC_LANG_PUSH(C++)
  # test compilation of simple program
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $CULA_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $CULA_CFLAGS"

  AC_TRY_LINK([#include <cula.hpp>],[culaStatus status = FakeculaInitialize();],
              have_cula=yes, have_cula=no)

  if test $have_cula = no; then
    AC_TRY_LINK([#include <cula_status.h>],[FakeculaInitialize();],
              have_cula=yes, have_cula=no)
  fi

  AC_MSG_RESULT($have_cula)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_cula" = xyes; then
    AC_DEFINE([HAVE_CULA], [1], [Define to 1 if you have the CULA library])
    [$1]
  else
    AC_MSG_WARN([Will compile without CULA GPU code])
    if test x"$CULA_LIB_PATH_64" = x; then
      AC_MSG_WARN([Please set the CULA_LIB_PATH_64/32 environment variable])
    fi
    CULA_CFLAGS=""
    CULA_LIBS=""
    [$2]
  fi

  AC_SUBST(CULA_CFLAGS)
  AC_SUBST(CULA_LIBS)
  AM_CONDITIONAL(HAVE_CULA, [test x"$have_cula" = xyes])

  AC_LANG_POP(C++)
])


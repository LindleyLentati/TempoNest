#
# SWIN_LIB_MLAPACK([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the high precision mlapack libraries.
#
# MLAPACK_CFLAGS - autoconfig variable with flags required for compiling
# MLAPACK_LIBS   - autoconfig variable with flags required for linking
# HAVE_MLAPACK   - automake conditional
# HAVE_MLAPACK   - pre-processor macro in config.h
#
# Some environment variables are required for the MLAPACK libraries to
# function, if they are not installed in standard locations
#
# $MLAPACK need to be pointing to the installation
# directory
#
# This macro tries to link a test program in the following way
#
#    -L$(MLAPACK)/lib  -lmblas_qd -lmlapack_qd -lmblas_dd -lmlapack_dd  -I$(MLAPACK)/include
#
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_MLAPACK],
[
  AC_PROVIDE([SWIN_LIB_MLAPACK])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])

  AC_MSG_CHECKING([for MLAPACK installation])

  MLAPACK_CFLAGS=""
  MLAPACK_LIBS=""

  MLAPACK_LIB=""

  if test x"$MLAPACK" != x; then
    MLAPACK_LIBS="-L$MLAPACK/lib -lmblas_qd -lmlapack_qd -lmblas_dd -lmlapack_dd"
  fi

  if test x"$MLAPACK" != x; then
    MLAPACK_CFLAGS="-I$MLAPACK/include"
  fi

  MLAPACK_LIBS="$MLAPACK_LIBS $MLAPACK_LIB"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $MLAPACK_LIBS"
  CFLAGS="$ac_save_CFLAGS $MLAPACK_CFLAGS"

  AC_LANG_PUSH(C++)
  # test compilation of simple program
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $MLAPACK_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $MLAPACK_CFLAGS"

  AC_TRY_LINK([#include "mpack/mblas_dd.h"],[Mlsame_dd("", "");],
              have_mlapack=yes, have_mlapack=no)


  AC_MSG_RESULT($have_mlapack)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_mlapack" = xyes; then
    AC_DEFINE([HAVE_MLAPACK], [1], [Define to 1 if you have the MLAPACK library])
    [$1]
  else
    AC_MSG_WARN([Will compile without MLAPACK code])
    if test x"$MLAPACK" = x; then
      AC_MSG_WARN([Please set the MLAPACK environment variable])
    fi
    MLAPACK_CFLAGS=""
    MLAPACK_LIBS=""
    [$2]
  fi

  AC_SUBST(MLAPACK_CFLAGS)
  AC_SUBST(MLAPACK_LIBS)
  AM_CONDITIONAL(HAVE_MLAPACK, [test x"$have_mlapack" = xyes])

  AC_LANG_POP(C++)
])


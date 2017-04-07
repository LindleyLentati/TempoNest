#
# SWIN_LIB_MULTINEST([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the package MultiNest
# by F. Feroz, M.P. hobson, and M. Bridges
#
# MULTINEST_CFLAGS - autoconfig variable with flags required for compiling
# MULTINEST_LIBS   - autoconfig variable with flags required for linking
# HAVE_MULTINEST   - automake conditional
# HAVE_MULTINEST   - pre-processor macro in config.h
#
# This macro tries to link a test program, by using
#
#    -L$MULTINEST_DIR -lnest3 -llapack -blas
#
# Notice that the environment variable MULTINEST_DIR is used. In the case the
# library 'libnest' is not installed in a default location, let MULTINEST_DIR
# point to the location of libnest
# MULTINEST_LIBS="$MULTINEST_LIBS -lchord -lnest3 -lellipsis"
#
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_MULTINEST],
[
  AC_PROVIDE([SWIN_LIB_MULTINEST])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])
  AC_REQUIRE([AX_PATH_GSL])

  AC_MSG_CHECKING([for MultiNest installation])

  MULTINEST_CFLAGS=""
  MULTINEST_LIBS=""

  if test x"$MULTINEST_DIR" != x; then
    MULTINEST_CFLAGS="-I$MULTINEST_DIR"
    MULTINEST_LIBS="-L$MULTINEST_DIR"
  fi

  MULTINEST_LIBS="$MULTINEST_LIBS -lchord -lnest3"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $MULTINEST_LIBS"
  CFLAGS="$ac_save_CFLAGS $MULTINEST_CFLAGS"

  AC_TRY_LINK([#include <stdio.h>],[__nested_MOD_nestrun(NULL);],
              have_multinest=yes, have_multinest=yes)

  if test $have_multinest = no; then
    AC_TRY_LINK([#include <stdio.h>],[nested_mp_nestrun_(NULL);],
              have_multinest=yes, have_multinest=no)
  fi

  AC_MSG_RESULT($have_multinest)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_multinest" = xyes; then
    AC_DEFINE([HAVE_MULTINEST], [1], [Define to 1 if you have the MULTINEST library])
    [$1]
  else
    AC_MSG_WARN([MULTINEST code will not be compiled])
    if test x"$MULTINEST_DIR" = x; then
      AC_MSG_WARN([Please set the MULTINEST_DIR environment variable])
    fi
    MULTINEST_CFLAGS=""
    MULTINEST_LIBS=""
    [$2]
  fi

  AC_SUBST(MULTINEST_CFLAGS)
  AC_SUBST(MULTINEST_LIBS)
  AM_CONDITIONAL(HAVE_MULTINEST, [test x"$have_multinest" = xyes])
])


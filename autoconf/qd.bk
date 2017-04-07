#
# SWIN_LIB_QDINSTALL([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the high precision qdinstall libraries.
#
# QDINSTALL_CFLAGS - autoconfig variable with flags required for compiling
# QDINSTALL_LIBS   - autoconfig variable with flags required for linking
# HAVE_QDINSTALL   - automake conditional
# HAVE_QDINSTALL   - pre-processor macro in config.h
#
# Some environment variables are required for the QDINSTALL libraries to
# function, if they are not installed in standard locations
#
# $QDINSTALL need to be pointing to the installation
# directory
#
# This macro tries to link a test program in the following way
#
#    -L$(QDINSTALL)/lib  -lmblas_qd -lqdinstall_qd -lmblas_dd -lqdinstall_dd  -I$(QDINSTALL)/include
#
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_QDINSTALL],
[
  AC_PROVIDE([SWIN_LIB_QDINSTALL])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])

  AC_MSG_CHECKING([for QDINSTALL installation])

  QDINSTALL_CFLAGS=""
  QDINSTALL_LIBS=""

  QDINSTALL_LIB=""

  if test x"$QDINSTALL" != x; then
    QDINSTALL_LIBS="-L$QDINSTALL/lib -lqd"
  fi

  if test x"$QDINSTALL" != x; then
    QDINSTALL_CFLAGS="-I$QDINSTALL/include"
  fi

  QDINSTALL_LIBS="$QDINSTALL_LIBS $QDINSTALL_LIB"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $QDINSTALL_LIBS"
  CFLAGS="$ac_save_CFLAGS $QDINSTALL_CFLAGS"

  AC_LANG_PUSH(C++)
  # test compilation of simple program
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $QDINSTALL_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $QDINSTALL_CFLAGS"

  AC_TRY_LINK([#include "qd/qd_real.h"],[qd_real rand();],
              have_qdinstall=yes, have_qdinstall=no)


  AC_MSG_RESULT($have_qdinstall)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_qdinstall" = xyes; then
    AC_DEFINE([HAVE_QDINSTALL], [1], [Define to 1 if you have the QDINSTALL library])
    [$1]
  else
    AC_MSG_WARN([Will compile without QDINSTALL code])
    if test x"$QDINSTALL" = x; then
      AC_MSG_WARN([Please set the QDINSTALL environment variable])
    fi
    QDINSTALL_CFLAGS=""
    QDINSTALL_LIBS=""
    [$2]
  fi

  AC_SUBST(QDINSTALL_CFLAGS)
  AC_SUBST(QDINSTALL_LIBS)
  AM_CONDITIONAL(HAVE_QDINSTALL, [test x"$have_qdinstall" = xyes])

  AC_LANG_POP(C++)
])


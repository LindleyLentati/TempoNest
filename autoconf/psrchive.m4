#
# SWIN_LIB_PSRCHIVE([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the psrchive libraries.
#
# PSRCHIVE_CFLAGS - autoconfig variable with flags required for compiling
# PSRCHIVE_LIBS   - autoconfig variable with flags required for linking
# HAVE_PSRCHIVE   - automake conditional
# HAVE_PSRCHIVE   - pre-processor macro in config.h
#
# Some environment variables are required for the PSRCHIVE libraries to
# function, if they are not installed in standard locations
#
# $PSRCHIVE need to be pointing to the installation
# directory
#
# This macro tries to link a test program in the following way
#
#    -L$(PSRCHIVE)/lib -lpsrbase -lpsrmore -lpsrutil  -I$(PSRCHIVE)/include
#
# TODO: Test whether the BLAS and LAPACK libraries really are required
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_PSRCHIVE],
[
  AC_PROVIDE([SWIN_LIB_PSRCHIVE])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])

  AC_MSG_CHECKING([for PSRCHIVE installation])

  PSRCHIVE_CFLAGS=""
  PSRCHIVE_LIBS=""

  PSRCHIVE_LIB="-lpsrbase -lpsrmore -lpsrutil -lstdc++  -lfftw3f -lcfitsio -lgfortran -ltempo2"

  if test x"$PSRCHIVE" != x; then
    PSRCHIVE_LIBS="-L$PSRCHIVE/lib"
  fi

  if test x"$PSRCHIVE" != x; then
    PSRCHIVE_CFLAGS="-I$PSRCHIVE/include"
  fi

  PSRCHIVE_LIBS="$PSRCHIVE_LIBS $PSRCHIVE_LIB"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $PSRCHIVE_LIBS"
  CFLAGS="$ac_save_CFLAGS $PSRCHIVE_CFLAGS"

  AC_LANG_PUSH(C++)
  # test compilation of simple program
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $PSRCHIVE_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $PSRCHIVE_CFLAGS"

  AC_TRY_LINK([#include "Pulsar/Archive.h"],[Pulsar::Archive::load("");],
              have_psrchive=yes, have_psrchive=yes)


  AC_MSG_RESULT($have_psrchive)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_psrchive" = xyes; then
    AC_DEFINE([HAVE_PSRCHIVE], [1], [Define to 1 if you have the PSRCHIVE library])
    [$1]
  else
    AC_MSG_WARN([Will compile without PSRCHIVE code])
    if test x"$PSRCHIVE" = x; then
      AC_MSG_WARN([Please set the PSRCHIVE environment variable])
    fi
    PSRCHIVE_CFLAGS=""
    PSRCHIVE_LIBS=""
    [$2]
  fi

  AC_SUBST(PSRCHIVE_CFLAGS)
  AC_SUBST(PSRCHIVE_LIBS)
  AM_CONDITIONAL(HAVE_PSRCHIVE, [test x"$have_psrchive" = xyes])

  AC_LANG_POP(C++)
])


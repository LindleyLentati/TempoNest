#
# SWIN_LIB_TEMPO2([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# This m4 macro checks availability of the tempo2 library by G.B. Hobbs and
# Russel Edwards
#
# TEMPO2_CFLAGS - autoconfig variable with flags required for compiling
# TEMPO2_LIBS   - autoconfig variable with flags required for linking
# HAVE_TEMPO2   - automake conditional
# HAVE_TEMPO2   - pre-processor macro in config.h
#
# This macro tries to link a test program, by using
#
#    -L$TEMPO2/lib -ltempo2
#
# Notice that the environment variable $TEMPO2 is used. This is standard for a
# Tempo2 installation.
#
# If the Tempo2 library is not installed in $TEMPO2/lib or a standard location,
# the variable $TEMPO2_LIB needs to be set.
#
# If the Tempo2 headers are not installed in $TEMPO2/include or a standard
# location, the variable $TEMPO2_INC needs to be set.
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_TEMPO2],
[
  AC_PROVIDE([SWIN_LIB_TEMPO2])

  AC_MSG_CHECKING([for TEMPO2 installation])

  TEMPO2_CFLAGS=""
  TEMPO2_LIBS=""

  TEMPO2_LIB=""

  if test x"$TEMPO2" != x; then
    TEMPO2_LIBS="-L$TEMPO2/lib -ltempo2"
  fi

  if test x"$TEMPO2" != x; then
    TEMPO2_CFLAGS="-I$TEMPO2/include"
  fi

  TEMPO2_LIBS="$TEMPO2_LIBS $TEMPO2_LIB"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $TEMPO2_LIBS"
  CFLAGS="$ac_save_CFLAGS $TEMPO2_CFLAGS"

  AC_LANG_PUSH(C++)
  # test compilation of simple program
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $TEMPO2_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $TEMPO2_CFLAGS"

  AC_TRY_LINK([#include "tempo2.h"],[readParfile("");],
              have_tempo2=yes, have_tempo2=yes)


  AC_MSG_RESULT($have_tempo2)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_tempo2" = xyes; then
    AC_DEFINE([HAVE_TEMPO2], [1], [Define to 1 if you have the TEMPO2 library])
    [$1]
  else
    AC_MSG_WARN([Will compile without TEMPO2 code])
    if test x"$TEMPO2" = x; then
      AC_MSG_WARN([Please set the TEMPO2 environment variable])
    fi
    TEMPO2_CFLAGS=""
    TEMPO2_LIBS=""
    [$2]
  fi

  AC_SUBST(TEMPO2_CFLAGS)
  AC_SUBST(TEMPO2_LIBS)
  AM_CONDITIONAL(HAVE_TEMPO2, [test x"$have_tempo2" = xyes])

  AC_LANG_POP(C++)
])


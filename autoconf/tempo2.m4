#
# SWIN_LIB_TEMPO2([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
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

  AC_MSG_CHECKING([for Tempo2 installation])

  TEMPO2_CFLAGS=""
  TEMPO2_LIBS=""

  if test x"$TEMPO2_INC" != x; then
    TEMPO2_CFLAGS="-I$TEMPO2_INC"
    TEMPO2_INCDIR="$TEMPO2_INC"
  else
    TEMPO2_CFLAGS="-I$TEMPO2/include"
    TEMPO2_INCDIR="$TEMPO2/include"
  fi

  if test x"$TEMPO2_LIB" != x; then
    TEMPO2_LIBS="-L$TEMPO2_LIB"
    TEMPO2_LIBDIR="$TEMPO2_LIB"
  else
    TEMPO2_LIBS="-L$TEMPO2/lib"
    TEMPO2_LIBDIR="$TEMPO2/lib"
  fi

  TEMPO2_LIBS="$TEMPO2_LIBS -ltempo2"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $TEMPO2_LIBS"
  CFLAGS="$ac_save_CFLAGS $TEMPO2_CFLAGS"

  AC_TRY_LINK([#include <stdio.h>],[readParfile("");],
              have_tempo2=yes, have_tempo2=yes)

  AC_MSG_RESULT($have_tempo2)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_tempo2" = xyes; then
    AC_DEFINE([HAVE_TEMPO2], [1], [Define to 1 if you have the TEMPO2 library])
    [$1]
  else
    AC_MSG_WARN([TEMPO2 code will not be compiled])
    if test x"$TEMPO2_LIB" = x; then
      AC_MSG_WARN([Please set the TEMPO2_LIB environment variable])
    fi
    if test x"$TEMPO2_INC" = x; then
      AC_MSG_WARN([Please set the TEMPO2_INC environment variable])
    fi
    TEMPO2_CFLAGS=""
    TEMPO2_LIBS=""
    TEMPO2_LIBDIR=""
    TEMPO2_INCDIR=""
    [$2]
  fi

  if test x"$have_tempo2" = xyes; then
    CFLAGS="$CFLAGS $TEMPO2_CFLAGS"
    LDFLAGS="$LDFLAGS $TEMPO2_LIBS"

    dnl In tempo2.h, the version is denoted as:
    dnl #define TEMPO2_h_VER "$Revision: 1.65 $"
    dnl This seems to be the only way to obtain the header version
    dnl Convert to float from the 11-th character of the string onwards
    dnl TODO: How to get the version number as a string here?
    AC_MSG_CHECKING([For for compatible tempo2.h])
    AC_RUN_IFELSE(
      [AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT
#include <tempo2.h>
#include <stdlib.h>
          ],
        [ #ifdef TEMPO2_h_VER
            if(atof(TEMPO2_h_VER+11) >= 1.65) {
                return 0;
            } else {
                return -1;
            }
          #else
            return -1;
          #endif
        ])
      ],
      [
        AC_MSG_RESULT([yes])
      ],
      [

        AC_RUN_IFELSE(
          [AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT
#include <tempo2.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
              ],
            [ #ifdef TEMPO2_h_VER
                if(strncmp(TEMPO2_h_VER, "$Revision$", 10) == 0) {
                    return 0;
                } else {
                    return 0;
                }
              #else
                return -1;
              #endif
            ])
          ],
          [
            AC_MSG_RESULT([cvs version])
            AC_MSG_WARN([[
========================================================================
CVS version of tempo2 detected. The tempo2 version cannot be reliably
determined from a CVS build. A tempo2 version >= 1.65 is required. If
the build does not succeed, please update tempo2.
========================================================================]])
          ],
          [

            AC_MSG_RESULT([not found or version too old (requires >= 1.65)])
            AC_MSG_WARN([Please install a more recent version of Tempo2])
            have_tempo2=no
            TEMPO2_CFLAGS=""
            TEMPO2_LIBS=""
            TEMPO2_LIBDIR=""
            TEMPO2_INCDIR=""
          ])
      ])
  fi

  AC_SUBST(TEMPO2_CFLAGS)
  AC_SUBST(TEMPO2_LIBS)
  AC_SUBST(TEMPO2_LIBDIR)
  AC_SUBST(TEMPO2_INCDIR)
  AM_CONDITIONAL(HAVE_TEMPO2, [test x"$have_tempo2" = xyes])

])


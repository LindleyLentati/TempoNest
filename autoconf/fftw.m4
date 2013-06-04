dnl @synopsis SWIN_LIB_FFTW
dnl 
AC_DEFUN([SWIN_LIB_FFTW],
[
  AC_PROVIDE([SWIN_LIB_FFTW])

  SWIN_PACKAGE_OPTIONS([fftw3])

  AC_MSG_CHECKING([for double-precision FFTW-3 library])

  if test "$have_fftw" != "user disabled"; then

    SWIN_PACKAGE_FIND([fftw3],[fftw3.h])
    SWIN_PACKAGE_TRY_COMPILE([fftw3],[#include <fftw3.h>])

    if test $have_fftw3 = yes; then
      SWIN_PACKAGE_FIND([fftw3],[libfftw3.*])
      SWIN_PACKAGE_TRY_LINK([fftw3],[#include <fftw3.h>],
                            [fftw_plan_dft_1d(0,0,0,FFTW_FORWARD,FFTW_ESTIMATE);],
                            [-lfftw3 -lm])
    fi

  else
    have_fftw3=no
  fi

  AC_MSG_RESULT([$have_fftw3])

  if test $have_fftw3 = yes; then
    AC_DEFINE(HAVE_FFTW3,1,[Define if the FFTW3 library is installed])
    FFTW_LIBS="$fftw3_LIBS $FFTW_LIBS"
    FFTW_CFLAGS="$fftw3_CFLAGS $FFTW_CFLAGS"
  fi

  AC_SUBST(FFTW_LIBS)
  AC_SUBST(FFTW_CFLAGS)

  AM_CONDITIONAL(HAVE_FFTW3,[test "$have_fftw3" = yes])

])


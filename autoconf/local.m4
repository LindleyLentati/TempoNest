#
# SWIN_LOCAL
#
# This m4 macro enables the user to optionally disable the use of packages
# installed in /usr/local or /sw and $PSRHOME/packages/$LOGIN_ARCH
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LOCAL],
[
  AC_PROVIDE([SWIN_LOCAL])

#
# Set up to use the standard package installation directories
#
AC_ARG_ENABLE([local],
              AC_HELP_STRING([--disable-local],
              [Don't use packages in /usr/local or /sw or /opt/local]))

if test x"$enable_local" != x"no"; then
  local_dirs="/usr/local /sw /opt/local"
  for dir in $local_dirs; do
    if test -d $dir; then
      CPPFLAGS="-I$dir/include $CPPFLAGS"
      LDFLAGS="-L$dir/lib $LDFLAGS"
    fi
  done
fi

#
# Set up to use the pulsar group package installation directory
#
AC_ARG_ENABLE([psrhome],
              AC_HELP_STRING([--disable-psrhome],
              [Don't use packages in $PSRHOME/packages/$LOGIN_ARCH]))

if test -d $PSRHOME/packages/$LOGIN_ARCH -a x"$enable_psrhome" != x"no"; then
  CPPFLAGS="-I$PSRHOME/packages/$LOGIN_ARCH/include $CPPFLAGS"
  LDFLAGS="-L$PSRHOME/packages/$LOGIN_ARCH/lib $LDFLAGS"
fi

if test -d $PSRHOME/$LOGIN_ARCH -a x"$enable_psrhome" != x"no"; then
  CPPFLAGS="-I$PSRHOME/$LOGIN_ARCH/include $CPPFLAGS"
  LDFLAGS="-L$PSRHOME/$LOGIN_ARCH/lib $LDFLAGS"
fi

])


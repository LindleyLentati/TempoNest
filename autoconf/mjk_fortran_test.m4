AC_DEFUN([MJK_FORTRAN_TEST],
[
  AC_PROVIDE([MJK_FORTRAN_TEST])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([AC_PROG_F77])
  AC_REQUIRE([AC_F77_WRAPPERS])

  AC_MSG_CHECKING([C and Fortran compilers are compatible])

  fortran_c_links=no
  AC_F77_FUNC(mktest)
  $ECHO "int main(int argc, char** argv){$mktest (); return 0;}" > conftest_c.c
  $ECHO '      SUBROUTINE mktest ()' > conftest_for.f
  $ECHO '      write(*,*) "TEST"' >> conftest_for.f
  $ECHO '      END' >> conftest_for.f

  $ECHO "$F77 $FFLAGS conftest_for.f -c -o conftest_for.o" >&5
  $F77 $FFLAGS conftest_for.f -c -o conftest_for.o >&5
  ex=$?
  if test $ex -eq 0  ; then
    $ECHO "$CC conftest_c.c conftest_for.o $CFLAGS $LDFLAGS $LIBS $FLIBS -o conftest_cfor" >&5
    $CC conftest_c.c conftest_for.o $CFLAGS $LDFLAGS $LIBS $FLIBS -o conftest_cfor >&5
    ex=$?
    if test $ex -eq 0  ; then
      fortran_c_links=yes
    else
      $ECHO "Failed program was:" >&5
      cat conftest_c.c >&5
      cat conftest_for.f >&5
    fi
  fi

  AC_MSG_RESULT($fortran_c_links)

  rm conftest_c.c >& /dev/null
  rm conftest_for.c >& /dev/null
  rm conftest_for.o >& /dev/null
  rm conftest_cfor >& /dev/null


])


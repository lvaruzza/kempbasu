dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([side], [0.5.0], [varuzza@gmail.com])

dnl AM_PROG_CC_C_O

AM_INIT_AUTOMAKE(side,0.5.0)

AC_HEADER_STDC

AC_PROG_CC

ACX_PTHREAD()
AC_CHECK_HEADERS(Judy.h)
AC_CHECK_LIB(Judy, JudyHSIns, ,echo "SIDE requires Judy to run. Source code provied in ext-lib directory"; exit 1)
PKG_CHECK_MODULES(GSL, gsl)
AC_SUBST(GSL_LIBS)
AC_SUBST(GSL_CFLAGS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

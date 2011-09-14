dnl Process this file with autoconf to produce a configure script.

AC_INIT([IMa2],[4.11.11])
AC_CONFIG_SRCDIR(src/ima_main.c)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign no-dependencies])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

dnl things required by automake
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_LANG(C)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AM_PROG_CC_C_O

dnl Check compiler features
AC_TYPE_SIZE_T
dnl AC_C_CONST
AC_C_VOLATILE
AC_C_INLINE
AC_C_CHAR_UNSIGNED

dnl AC_FUNC_ALLOCA
AC_FUNC_VPRINTF

dnl
dnl IMa2 configuration
dnl

case "$host" in
    *-*-linux*) 
        platform_interface=linux
        ;;
    *-*-*darwin*) 
        platform_interface=darwin
        ;;
esac
 
if test "$platform_interface" = "linux" ; then
  AC_DEFINE([LINUX], 1, [Define for linux.])
fi

if test "$platform_interface" = "darwin" ; then
  AC_DEFINE([DARWIN], 1, [Define for mac.])
fi

AS_IF([test "x$STRUCTURAMA" = x], 
      [AC_DEFINE_UNQUOTED([STRUCTURAMA],["/usr/local/bin"],[Path of Structurama])])

dnl
dnl IMa2 configuration
dnl

dnl
AC_CONFIG_FILES([src/Makefile Makefile])
AC_OUTPUT




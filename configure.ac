# Process this file with autoconf to produce a configure script.
AC_INIT([sig], [1.0], [edeveaud@pasteur.fr])
AC_CONFIG_SRCDIR([src/main.c])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([src/config.h])

# Checks for programs.
AC_PROG_CC
AC_CHECK_PROG(POD2MAN, pod2man, pod2man, :)
AC_CHECK_PROG(PODCHECKER, podchecker, podchecker, :)

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([errno.h stdlib.h string.h unistd.h ctype.h sys/types.h])
AC_CHECK_HEADERS([regex.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([isascii regcomp strerror strrchr])

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile test/Makefile])
AC_OUTPUT

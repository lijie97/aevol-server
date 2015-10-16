#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Require version 2.65 (or newer) of the autotools
AC_PREREQ([2.65])

# Provide FULL-PACKAGE-NAME, VERSION and BUG-REPORT-ADDRESS
AC_INIT([AEVOL], [5.0-beta6], [aevol-bugs@lists.gforge.inria.fr], , [www.aevol.fr])

# Set C++ as the main language
AC_LANG([C++])

# Change default value for CXXFLAGS
test_CXXFLAGS=${CXXFLAGS+set}
if test "$test_CXXFLAGS" != set; then
  CXXFLAGS=""
fi

# Specify the location of additional local Autoconf macros
AC_CONFIG_MACRO_DIR([m4])

# We want to use automake
AM_INIT_AUTOMAKE

# Detect C++11 or C++14
AX_CXX_COMPILE_STDCXX_14([noext])
if test -z "$HAVE_CXX14" ; then
   AX_CXX_COMPILE_STDCXX_11([noext])
fi

#if test -z "$HAVE_CXX14";
#fi

# Check the existance of an arbitrary source file
AC_CONFIG_SRCDIR([src/aevol_run.cpp])

# Build config.h using config.h.in as a template
#AC_CONFIG_HEADERS([config.h:config.h.in])

# Checks for programs
AC_PROG_CXX([clang-omp++ clang++ g++])
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_OPENMP
AM_PROG_CC_C_O

# Print warning if openmp not found
AC_CHECK_HEADER([omp.h], , [AC_MSG_WARN([openmp not found - continuing without openmp support])])

# Initialize configure-time defined flags
AC_SUBST(AEVOLCPPFLAGS, "") # C PreProcessor flags
AC_SUBST(AEVOLCXXFLAGS, "$OPENMP_CXXFLAGS") # C++ Compile flags
AC_SUBST(AEVOLLDFLAGS, "$OPENMP_CXXFLAGS")  # Linker flags

# Checks for libraries.
AC_CHECK_LIB([m], [cos], , [AC_MSG_ERROR([PACKAGE_NAME requires the math library, please install it.])])
AC_CHECK_LIB([z], [gzread], , [AC_MSG_ERROR([PACKAGE_NAME requires the zlib development library, please install it.])])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h limits.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([floor memset mkdir pow rint sqrt gzread])

# ****************************************************************************
#                                  OPTIONS
#
#   Supported options [default]:
#     --with-x                          [yes]
#     --enable-optim                    [enabled]
#     --enable-normalized-fitness       [disabled]
#     --enable-mtperiod=period          [disabled]
#     --enable-trivialjumps=jumpsize    [disabled]
#     --enable-devel                    [disabled]
#     --enable-debug                    [disabled]
#     --enable-profiling                [disabled]
#
# ****************************************************************************

#   ************************ --without-x option *********************
#
#   When this option is set, typing make will produce aevol and
#   most post-treatment programs, but not aevol_X11.
#   This is useful if we want to run computations on a cluster where the
#   X libraries are not present. If, on the contrary, we type ./configure
#   without this option, then typing make will produce aevol_X11 and
#   view_generation_X11 (and the other post-treatment programs).

AC_MSG_CHECKING([whether to enable graphical outputs])
AC_ARG_WITH(x, [  --without-x to disable graphical output], , with_x=yes)
AC_MSG_RESULT($with_x)
AM_CONDITIONAL(WITH_X, test x"$with_x" != xno)
if test "x$with_x" != xno ; then
   AC_DEFINE(__X11)
   AC_PATH_XTRA
   # AC_PATH_XTRA updates X_LIBS and X_CFLAGS but not CPPFLAGS, CXXFLAGS or LDFLAGS.
   # Thus we need to add X_CFLAGS to AEVOLCPPFLAGS, and X_LIBS to AEVOLLDFLAGS.
   # Moreover, on MacOS, AC_CHECK_LIB needs correct CPPFLAGS to find X11. But CPPFLAGS is normally a user
   # variable, which we should not change. Thus we save its previous value to restore it after AC_CHECK_LIB.
   ac_cppflags_save="$CPPFLAGS"
   CPPFLAGS+=" $X_CFLAGS"
   AEVOLCPPFLAGS+=" $X_CFLAGS"
   AEVOLLDFLAGS+=" $X_LIBS"
   AC_CHECK_LIB(X11, XOpenDisplay, , [AC_MSG_ERROR([the X11 library could not be found on your system. Please install it or re-run with the --without-x option])],$X_LIBS)
   AC_CHECK_HEADERS([X11/Xlib.h X11/Xutil.h X11/keysym.h X11/Xatom.h X11/cursorfont.h])
   CPPFLAGS="$ac_cppflags_save"
else
   AC_DEFINE(__NO_X)
fi

#  ******************* --enable-optim option ***************************
#
#  With this option, the programs are compiled with the
#  optimization options (-O3 - fexpensive-optimizations...).
#

AC_MSG_CHECKING([whether to use compile time optimization])
AC_ARG_ENABLE(optim, [  --disable-optim turn off compile time optimization], , enable_optim=yes)
AC_MSG_RESULT($enable_optim)
AM_CONDITIONAL(WITH_OPTIM, test x"$enable_optim" != xno)
if test "x$enable_optim" = "xyes" ; then
  AEVOLCXXFLAGS+=" -O3 -ffast-math -pipe -fexpensive-optimizations"
else
  AEVOLCXXFLAGS+=" -O0"
fi

#  ******************* --enable-normalized-fitness option ***************************
#
#  With this option, the NORMALIZED_FITNESS flag is defined,
#  allowing a different fitness calculation

AC_MSG_CHECKING([whether to use normalized fitness])
AC_ARG_ENABLE(normalized-fitness, [  --enable-normalized-fitness to enable fitness normalization], , enable_normalized_fitness=no)
AC_MSG_RESULT($enable_normalized_fitness)
AM_CONDITIONAL(WITH_NORMALIZED_FITNESS, test x"$enable_normalized_fitness" != xno)
if test "x$enable_normalized_fitness" = "xyes" ; then
   AC_DEFINE(NORMALIZED_FITNESS)
fi

#  ************************* --enable-mtperiod=period option *******************************
#
#  This option sets a particular period for the Mersenne Twister.
#  Default is --enable-mtperiod=607 corresponding to a period of 2^607.
#  Possible values are 607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049 and 216091
#

AC_MSG_CHECKING([whether to use a particular period for the Mersenne Twister])

AC_ARG_ENABLE(mtperiod,
[  --enable-mtperiod=period
        Set the Mersenne Twister period to 2^period (default is 2^607).

        Possible values are 607, 1279, 2281, 4253, 11213, 19937, 44497,
          86243, 132049 and 216091

        Consider using in conjunction with --enable-trivialjumps for
          large periods
],
period=$enableval,
period=607)

if  test $period != 607 && test $period != 1279 && test $period != 2281 && test $period != 4253 && test $period != 11213 && test $period != 19937 &&
    test $period != 44497 && test $period != 86243 && test $period != 132049 && test $period != 216091 ; then
  AC_MSG_RESULT([FAIL])
  AC_MSG_ERROR([$period is not a valid Mersenne Twister period])
fi

AC_MSG_RESULT($period)
AC_DEFINE_UNQUOTED([SFMT_MEXP], [$period])

#  ************************* --enable-trivialjumps=jumpsize option *******************************
#
#  When this option is set, a trivial jump algorithm will be used instead of the polynomial-based method.
#  This is recommended for large Mersenne Twister periods when huge jumps are not needed since
#  the polynomial-based method takes more time in that case.
#
#  For the standard MT (period 2^19937), it takes approx. the same time to draw a few million numbers
#  than it does to make the same jump using the polynomial-based method.
#
#  The jumpsize argument allows to set the desired jump size. The default value is 1000 which is enough to
#  avoid overlaps in most cases.
#

AC_MSG_CHECKING([whether to use trivial Mersenne Twister jumps])
AC_ARG_ENABLE(trivialjumps,
[  --enable-trivialjumps=jumpsize
        When this option is set, a trivial jump algorithm will be used instead
        of the polynomial-based method.

        This is recommended for large Mersenne Twister periods when huge jumps
        are not needed since the polynomial-based method takes more time in that
        case.

        For the standard MT (period 2^19937), it takes approx. the same time to
        draw a few million numbers than it does to make the same jump using the
        polynomial-based method.

        The jumpsize argument allows to set the desired jump size. The default
        value is 1000 which is enough to avoid overlaps in most cases.
],
trivialjumps=$enableval,
trivialjumps=no)
if test "x$trivialjumps" != "xyes" && test "x$trivialjumps" != "xno" ; then
  AC_MSG_RESULT([yes ($trivialjumps)])
else
  AC_MSG_RESULT([$trivialjumps])
fi
AM_CONDITIONAL(ENABLE_FORCEPJ, test x"$trivialjumps" != xno)
if test "x$trivialjumps" = "xyes" ; then
  AC_DEFINE([TRIVIAL_METHOD_JUMP_SIZE], [1000])
elif test "x$trivialjumps" != "xno" ; then
  AC_DEFINE_UNQUOTED([TRIVIAL_METHOD_JUMP_SIZE], [$trivialjumps])
fi


#  ******************* --enable-devel option ***************************
#
#  This enables developers' features:
#    adds a call to make clean at the end of the configure script.
#    enables doxygen if available
#

AC_MSG_CHECKING([whether to use devel features])
AC_ARG_ENABLE(devel, [  --enable-devel turn on development features], , enable_devel=no)
AC_MSG_RESULT($enable_devel)
AM_CONDITIONAL(ENABLE_DEVEL, test x"$enable_devel" != xno)
if test "x$enable_devel" = "xyes" ; then
  # Make clean at the end of the configuration process
  # (usefull when several configurations are done consecutively, typically during development)
  AC_CONFIG_COMMANDS(make, [make clean])

  # Check for Doxygen
  AC_CHECK_PROGS([HAVE_DOXYGEN], [doxygen])
  if test -z "$HAVE_DOXYGEN";
     then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
  fi
else
  # Output a warning message about multi-configuration and make clean
  AC_CONFIG_COMMANDS([user_warning], [AC_MSG_NOTICE([Consider using make clean if you are re-configuring.])])
fi
AM_CONDITIONAL([WANT_DOXYGEN], [test -n "$HAVE_DOXYGEN"])


#  ******************* --enable-debug option ***************************
#
#  With this option, the DEBUG flag is defined, and the programs are
#  compiled with debugging options (-g -Wall).
#

AC_MSG_CHECKING([whether to use debugging features])
AC_ARG_ENABLE(debug, [  --enable-debug to enable degugging features], , enable_debug=no)
AC_MSG_RESULT($enable_debug)
AM_CONDITIONAL(WITH_DEBUG, test x"$enable_debug" != xno)
if test "x$enable_debug" = "xyes" ; then
   AC_DEFINE(DEBUG)
   AEVOLCXXFLAGS+=" -g -Wall"
else
   AC_DEFINE(NDEBUG)
fi


#  ******************* --enable-profiling option ***************************
#
#  With this option, the -pg flag is defined.
#

AC_MSG_CHECKING([whether to enable profiling])
AC_ARG_ENABLE(profiling, [  --enable-profiling to enable profiling], , enable_profiling=no)
AC_MSG_RESULT($enable_profiling)
AM_CONDITIONAL(WITH_PROFILING, test x"$enable_profiling" != xno)
if test "x$enable_profiling" = "xyes" ; then
   AEVOLCXXFLAGS+=" -pg"
   AEVOLLDFLAGS+=" -pg"
fi

#  ************************* back to GENERAL stuff *******************************

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/libaevol/Makefile
                 src/libaevol/SFMT-src-1.4/Makefile
                 src/post_treatments/Makefile
                 doc/Makefile])
AC_OUTPUT

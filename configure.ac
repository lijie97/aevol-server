#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Require version 2.65 (or newer) of the autotools
AC_PREREQ([2.65])

# Provide FULL-PACKAGE-NAME, VERSION and BUG-REPORT-ADDRESS
AC_INIT([AEVOL], [4.0.0.0], [aevol-bug-report@lists.gforge.liris.cnrs.fr], , [www.aevol.fr])

# We want to use automake
AM_INIT_AUTOMAKE

# Check the existance of an arbitrary source file
AC_CONFIG_SRCDIR([src/aevol.cpp])

# Build config.h using config.h.in as a template
#AC_CONFIG_HEADERS([config.h:config.h.in])

# Checks for programs (c anc c++ compilers)
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_CC_C_O
#AC_PROG_DOXYGEN

# Checks for libraries.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([z], [gzread])

# Checks for header files.
AC_PATH_X
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

# Make clean at the end of the configuration process
# (usefull when several configurations are done consecutively, typically during development)
AC_CONFIG_COMMANDS(make, [make clean])








dnl ****************************************************************************
dnl                                  OPTIONS
dnl
dnl             Default : without debug (with optimizations)
dnl                       with    x
dnl                       without regulation
dnl                       with    in2p3
dnl ****************************************************************************



dnl  ******************* check option consistency ***************************
AC_MSG_CHECKING([checking parameter consistency])
if test "x$with_in2p3" = "xyes" ; then
  if test "x$with_debug" = "xyes" ; then
    AC_MSG_RESULT([FAIL])
    AC_MSG_ERROR([debug and in2p3 are not compatible])
  fi
  if test "x$with_x" != xno ; then
    AC_MSG_RESULT([FAIL])
    AC_MSG_ERROR([x and in2p3 are not compatible])
  fi
fi
AC_MSG_RESULT([OK])



dnl  ******************* --with-debug option ***************************
dnl
dnl  With this option, the DEBUG flag is defined, and the programs are 
dnl  compiled with debugging options (-g -Wall). Without this option, the 
dnl  NDEBUG flag is defined, and the programs are compiled with the
dnl  optimization options (-O3 - fexpensive-optimizations...). 


AC_MSG_CHECKING([whether to use debugging features])
AC_ARG_WITH(debug, [  --with-debug to enable degugging features], , with_debug=no)
AC_MSG_RESULT($with_debug)
AM_CONDITIONAL(WITH_DEBUG, test x"$with_debug" != xno)
if test "x$with_debug" = "xyes" ; then
   AC_DEFINE(DEBUG)
   CXXFLAGS="-g -Wall"
else
   AC_DEFINE(NDEBUG)
   CXXFLAGS="-O3 -ffast-math -pipe -fexpensive-optimizations"
fi






dnl   ************************ --without-x option *********************
dnl
dnl   When this option is set, typing make will produce aevol, raevol and
dnl   the post-treatment programs, but not aevol_X11 or raevol_X11. 
dnl   This is useful if we want to run computations on a cluster where the 
dnl   X libraries are not present. If, on the contrary, we type ./configure 
dnl   without this option, then typing make will produce aevol_X11, 
dnl   raevol_X11 and view_generation_X11 (and the other post-treatments 
dnl   programs).

AC_MSG_CHECKING([whether to enable graphical outputs])
AC_ARG_WITH(x, [  --without-x to disable graphical output], , with_x=yes)
AC_MSG_RESULT($with_x)
AM_CONDITIONAL(WITH_X, test x"$with_x" != xno)
if test "x$with_x" != xno ; then
   AC_DEFINE(__X11)
   AC_PATH_XTRA
   LDFLAGS=($X_LIBS)
   AC_CHECK_LIB(X11, XOpenDisplay, , ,$X_LIBS)
   AC_CHECK_LIB(Xi, XGetExtensionVersion, , ,$X_LIBS)
   AC_CHECK_LIB(Xmu, XmuLookupStandardColormap, , ,$X_LIBS)
   AC_CHECK_HEADERS([X11/Xlib.h X11/Xutil.h X11/keysym.h X11/Xatom.h X11/cursorfont.h])
else
   AC_DEFINE(__NO_X)
fi

dnl  ******************* --with-normalized-fitness option ***************************
dnl
dnl  With this option, the NORMALIZED_FITNESS flag is defined, 
dnl  allowing a different fitness calculation 

AC_MSG_CHECKING([whether to use normalized fitness])
AC_ARG_WITH(normalized-fitness, [  --with-normalized-fitness to enable fitness normalization], , with_normalized_fitness=no)
AC_MSG_RESULT($with_normalized_fitness)
AM_CONDITIONAL(WITH_NORMALIZED_FITNESS, test x"$with_normalized_fitness" != xno)
if test "x$with_normalized_fitness" = "xyes" ; then
   AC_DEFINE(NORMALIZED_FITNESS)
fi



dnl  ************************* --enable-mtperiod=period option *******************************
dnl
dnl  This option sets a particular period for the Mersenne Twister.
dnl  Default is --enable-mtperiod=607 corresponding to a period of 2^607.
dnl  Possible values are 607, 1279, 2281, 4253, 11213, 19937, 44497, 86243, 132049 and 216091
dnl

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



dnl  ************************* --enable-trivialjumps=jumpsize option *******************************
dnl
dnl  When this option is set, a trivial jump algorithm will be used instead of the polynomial-based method.
dnl  This is recommended for large Mersenne Twister periods when huge jumps are not needed since
dnl  the polynomial-based method takes more time in that case.
dnl  
dnl  For the standard MT (period 2^19937), it takes approx. the same time to draw a few million numbers
dnl  than it does to make the same jump using the polynomial-based method.
dnl
dnl  The jumpsize argument allows to set the desired jump size. The default value is 1000 which is enough to
dnl  avoid overlaps in most cases.
dnl

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



dnl  ************************* --with-raevol option *******************************
dnl
dnl  With this option, regulation is enabled (RAevol) and the __REGUL flag is defined
dnl

AC_MSG_CHECKING([whether to use raevol])
AC_ARG_WITH(raevol, [  --with-raevol to enable regulation], raevol=$withval, raevol=no)
AC_MSG_RESULT($raevol)
AM_CONDITIONAL(WITH_RAEVOL, test x"$raevol" != xno)
if test "x$raevol" = "xyes" ; then
   AC_DEFINE(__REGUL)
fi




dnl  ******************* --with-in2p3 option ***************************
dnl
dnl  With this option, the IN2P3 flag is defined

AC_MSG_CHECKING([whether to use IN2P3 specific features])
AC_ARG_WITH(in2p3, [  --with-in2p3 to enable IN2P3 specific features], , with_in2p3=no)
AC_MSG_RESULT($with_in2p3)
AM_CONDITIONAL(WITH_IN2P3, test x"$with_in2p3" != xno)
if test "x$with_in2p3" = "xyes" ; then
   AC_DEFINE(__IN2P3)
fi




dnl  ************************* back to GENERAL stuff *******************************

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/libaevol/Makefile
                 src/libaevol/SFMT-src-1.4/Makefile
                 src/post_treatments/Makefile])
AC_OUTPUT

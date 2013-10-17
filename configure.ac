#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Require version 2.65 (or newer) of the autotools
AC_PREREQ([2.65])

# Provide FULL-PACKAGE-NAME, VERSION and BUG-REPORT-ADDRESS
AC_INIT([AEVOL], [4.2], [aevol-bug-report@lists.gforge.liris.cnrs.fr], , [www.aevol.fr])

# We want to use automake
AM_INIT_AUTOMAKE

# Check the existance of an arbitrary source file
AC_CONFIG_SRCDIR([src/aevol_run.cpp])

# Build config.h using config.h.in as a template
#AC_CONFIG_HEADERS([config.h:config.h.in])

# Checks for programs (c anc c++ compilers)
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_PROG_CC_C_O
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])  




# Checks for libraries.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([z], [gzread])


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
#             Default : without regulation
#                       without debug (with optimizations)
#                       with    x
#                       with    cppunit (requires cppunit)
#                       without normalized fitness
#                       Mersenne Twister (MT) period of 2^607.
#                       MT jumps done through polynomial-based method
#                       with    in2p3
# ****************************************************************************



#  ******************* check option consistency ***************************
AC_MSG_CHECKING([checking parameter consistency])
if test "x$enable_in2p3" = "xyes" ; then
  if test "x$enable_debug" = "xyes" ; then
    AC_MSG_RESULT([FAIL])
    AC_MSG_ERROR([debug and in2p3 are not compatible])
  fi
  if test "x$with_x" != xno ; then
    AC_MSG_RESULT([FAIL])
    AC_MSG_ERROR([x and in2p3 are not compatible])
  fi
fi
AC_MSG_RESULT([OK])



#  ******************* --enable-debug option ***************************
#
#  With this option, the DEBUG flag is defined, and the programs are 
#  compiled with debugging options (-g -Wall). Without this option, the 
#  NDEBUG flag is defined, and the programs are compiled with the
#  optimization options (-O3 - fexpensive-optimizations...).
#

AC_MSG_CHECKING([whether to use debugging features])
AC_ARG_ENABLE(debug, [  --enable-debug to enable degugging features], , enable_debug=no)
AC_MSG_RESULT($enable_debug)
AM_CONDITIONAL(WITH_DEBUG, test x"$enable_debug" != xno)
if test "x$enable_debug" = "xyes" ; then
   AC_DEFINE(DEBUG)
   CXXFLAGS="-g -Wall"
   CPPFLAGS="-g -Wall"
else
   AC_DEFINE(NDEBUG)
   CXXFLAGS="-O3 -ffast-math -pipe -fexpensive-optimizations"
   CPPFLAGS="-O3 -ffast-math -pipe -fexpensive-optimizations"
fi



#  ******************* --enable-devel option ***************************
#
#  This enables developers' features. For now it only adds a call to
#  make clean at the end of the configure script.
#

AC_MSG_CHECKING([whether to use devel features])
AC_ARG_ENABLE(devel, [  --enable-devel], , enable_devel=no)
AC_MSG_RESULT($enable_devel)
AM_CONDITIONAL(ENABLE_DEVEL, test x"$enable_devel" != xno)
if test "x$enable_devel" = "xyes" ; then
  # Make clean at the end of the configuration process
  # (usefull when several configurations are done consecutively, typically during development)
  AC_CONFIG_COMMANDS(make, [make clean])
fi


#  **************** --without-cppunit option *************************
#
#  When this option is set, the test library will not be built.
#  This is useful if cppunit was not installed. However, developers
#  should install cppunit, configure with cppunit and run make check
#  to perform unit testing.

AC_MSG_CHECKING([whether to build the test library (requires cppunit)])
AC_ARG_WITH(cppunit, [  --without-cppunit to disable unit testing], , with_cppunit=yes)
AC_MSG_RESULT($with_cppunit)
AM_CONDITIONAL(WITH_CPPUNIT, test x"$with_cppunit" != xno)
if test "x$with_cppunit" != xno ; then
   cppunit_missing=no
   AM_PATH_CPPUNIT(1.10.2, , cppunit_missing=yes)
   if test "x$cppunit_missing" != xno ; then
      AC_MSG_ERROR([cppunit not found, make will fail - please install cppunit, or run ./configure --without-cppunit]) 
   fi 
   # AC_CHECK_LIB does not work with C++ libraries, this is why we use AM_PATH_CPPUNIT
   # This macro will check whether a Cppunit version >= 1.10.2 is installed
   # and, if, yes, will define CPPUNIT_CFLAGS and CPPUNIT_LIBS.
   # This macro is defined in the cppunit.m4 file provided by the cppunit package
   # when it is installed. Hence aclocal must be run before autoconf to have the
   # macro known by autoconf.  If cppunit is installed in a non-standard
   # directory, you need to specify where to find cppunit.m4. For example, if
   # cppunit is installed under /usr/local, run:
   #   ACLOCAL="aclocal -I /usr/local/share/aclocal" autoreconf -if
fi





#   ************************ --without-x option *********************
#
#   When this option is set, typing make will produce aevol, raevol and
#   the post-treatment programs, but not aevol_X11 or raevol_X11. 
#   This is useful if we want to run computations on a cluster where the 
#   X libraries are not present. If, on the contrary, we type ./configure 
#   without this option, then typing make will produce aevol_X11, 
#   raevol_X11 and view_generation_X11 (and the other post-treatments 
#   programs).

AC_MSG_CHECKING([whether to enable graphical outputs])
AC_ARG_WITH(x, [  --without-x to disable graphical output], , with_x=yes)
AC_MSG_RESULT($with_x)
AM_CONDITIONAL(WITH_X, test x"$with_x" != xno)
if test "x$with_x" != xno ; then
   AC_DEFINE(__X11)
   AC_PATH_XTRA
   CXXFLAGS="$CXXFLAGS $X_CFLAGS"
   CPPFLAGS="$CPPFLAGS $X_CFLAGS"
   LDFLAGS="$LDFLAGS $X_LIBS"
   AC_CHECK_LIB(X11, XOpenDisplay, , ,$X_LIBS)
   AC_CHECK_LIB(Xi, XGetExtensionVersion, , ,$X_LIBS)
   AC_CHECK_LIB(Xmu, XmuLookupStandardColormap, , ,$X_LIBS)
   AC_CHECK_HEADERS([X11/Xlib.h X11/Xutil.h X11/keysym.h X11/Xatom.h X11/cursorfont.h])
else
   AC_DEFINE(__NO_X)
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



#  ************************* --enable-raevol option *******************************
#
#  With this option, regulation is enabled (RAevol) and the __REGUL flag is defined
#

AC_MSG_CHECKING([whether to use raevol])
AC_ARG_ENABLE(raevol, [  --enable-raevol to enable regulation], raevol=$enableval, raevol=no)
AC_MSG_RESULT($raevol)
AM_CONDITIONAL(WITH_RAEVOL, test x"$raevol" != xno)
if test "x$raevol" = "xyes" ; then
   AC_DEFINE(__REGUL)
fi




#  ******************* --enable-in2p3 option ***************************
#
#  With this option, the IN2P3 flag is defined

AC_MSG_CHECKING([whether to use IN2P3 specific features])
AC_ARG_ENABLE(in2p3, [  --enable-in2p3 to enable IN2P3 specific features], , enable_in2p3=no)
AC_MSG_RESULT($enable_in2p3)
AM_CONDITIONAL(WITH_IN2P3, test x"$enable_in2p3" != xno)
if test "x$enable_in2p3" = "xyes" ; then
   AC_DEFINE(__IN2P3)
fi




#  ************************* back to GENERAL stuff *******************************

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/libaevol/Makefile
                 src/libaevol/SFMT-src-1.4/Makefile
                 src/libaevol/tests/Makefile
                 src/post_treatments/Makefile
                 tests/Makefile
                 doc/Makefile])
AC_OUTPUT

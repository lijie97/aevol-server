#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Require version 2.65 (or newer) of the autotools
AC_PREREQ([2.65])

# Provide FULL-PACKAGE-NAME, VERSION and BUG-REPORT-ADDRESS
AC_INIT([AEVOL], [3.1.0], aevol-bug-report@lists.gforge.liris.cnrs.fr)

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
AC_CHECK_FUNCS([floor memset mkdir pow rint sqrt])








dnl ****************************************************************************
dnl                                  OPTIONS
dnl
dnl             Default : without debug (with optimizations)
dnl                       with xlib
dnl                       without regulation
dnl                       with in2p3
dnl ****************************************************************************



dnl  ******************* check option consistency ***************************
AC_MSG_CHECKING([checking parameter consistency])
if test "x$with_in2p3" = "xyes" ; then
if test "x$with_debug" = "xyes" ; then
   AC_ERROR([debug and in2p3 are not compatible])
fi
if test "x$with_xlib" != xno ; then
   AC_ERROR([xlib and in2p3 are not compatible])
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
AC_ARG_WITH(debug, [--with-debug to enable degugging features],
		 with_debug=$withval,
		 with_debug=no)
AC_MSG_RESULT($with_debug)
AM_CONDITIONAL(WITH_DEBUG, test x"$with_debug" != xno)
if test "x$with_debug" = "xyes" ; then
   AC_DEFINE(DEBUG)
   CXXFLAGS="-g -Wall"
else
   AC_DEFINE(NDEBUG)
   CXXFLAGS="-O3 -ffast-math -pipe -fexpensive-optimizations"
fi





dnl   ************************ --without-xlib option *********************
dnl
dnl   When this option is set, typing make will produce aevol, raevol and
dnl   the post-treatment programs, but not aevol_X11 or raevol_X11. 
dnl   This is useful if we want to run computations on a cluster where the 
dnl   X libraries are not present. If, on the contrary, we type ./configure 
dnl   without this option, then typing make will produce aevol_X11, 
dnl   raevol_X11 and view_generation_X11 (and the other post-treatments 
dnl   programs).



AC_MSG_CHECKING([whether to enable graphical outputs])
AC_ARG_WITH(xlib, [--without-xlib to disable xlib usage completely],
		 with_xlib=$withval,
		 with_xlib=yes)
AC_MSG_RESULT($with_xlib)
AM_CONDITIONAL(WITH_X, test x"$with_xlib" != xno)
if test "x$with_xlib" != xno ; then
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



dnl  ************************* --with-raevol option *******************************
dnl
dnl  With this option, regulation is enabled (RAevol) and the __REGUL flag is defined
dnl


AC_MSG_CHECKING([whether to use raevol])
AC_ARG_WITH(raevol, [--with-raevol to enable regulation], raevol=$withval, raevol=no)
AC_MSG_RESULT($raevol)
AM_CONDITIONAL(WITH_RAEVOL, test x"$raevol" != xno)
if test "x$raevol" = "xyes" ; then
   AC_DEFINE(__REGUL)
fi




dnl  ******************* --with-in2p3 option ***************************
dnl
dnl  With this option, the IN2P3 flag is defined


AC_MSG_CHECKING([whether to use IN2P3 specific features])
AC_ARG_WITH(in2p3, [--with-in2p3 to enable IN2P3 specific features],
		 with_in2p3=$withval,
		 with_in2p3=no)
AC_MSG_RESULT($with_in2p3)
AM_CONDITIONAL(WITH_IN2P3, test x"$with_in2p3" != xno)
if test "x$with_in2p3" = "xyes" ; then
   AC_DEFINE(__IN2P3)
fi




dnl  ************************* back to GENERAL stuff *******************************

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/lib-aevol/Makefile
                 src/post_treatments/Makefile])
AC_OUTPUT

############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)

SUBDIRS = libaevol post_treatments

if WITH_GTEST
  SUBDIRS += tests
endif # WITH_GTEST

AM_CPPFLAGS += -I$(top_srcdir)/src/libaevol

############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS


####################################
# Add hardening option to compiler #
####################################
# This replaces insecure unlimited length buffer function calls with
# length-limited ones (required for debian and not a bad idea anyway).
AM_CPPFLAGS += -D_FORTIFY_SOURCE=2




############################################
# Set programs to build and library to use #
############################################


##### The solution below is no longer used since it generates a dependency problem #####
#AM_LDFLAGS = -L$(top_srcdir)/src/libaevol
#  LDADD = -laevol
#####  #####

LDADD = $(top_srcdir)/src/libaevol/libaevol.a $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a


bin_PROGRAMS = aevol_run aevol_create aevol_modify aevol_propagate
CLEANFILES = aevol_run aevol_create aevol_modify aevol_propagate

aevol_propagate_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS)
aevol_propagate_LDFLAGS = $(AM_LDFLAGS) $(BOOST_LDFLAGS)
aevol_propagate_LDADD = $(LDADD) $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)

aevol_run_SOURCES = aevol_run.cpp
aevol_create_SOURCES = aevol_create.cpp
aevol_modify_SOURCES = aevol_modify.cpp
aevol_propagate_SOURCES = aevol_propagate.cpp

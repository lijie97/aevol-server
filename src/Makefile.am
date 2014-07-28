############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)

SUBDIRS = libaevol post_treatments

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



if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = aevol_run_X11_debug
else # WITHOUT_DEBUG
  bin_PROGRAMS = aevol_run_X11
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = aevol_run_debug
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = aevol_run_IN2P3
else # WITHOUT_IN2P3
  bin_PROGRAMS = aevol_run
endif # IN2P3
endif # DEBUG
endif # X

bin_PROGRAMS += aevol_create aevol_modify aevol_propagate

CLEANFILES = aevol_run aevol_run_IN2P3 aevol_run_debug aevol_run_X11 aevol_run_X11_debug
CLEANFILES += aevol_create aevol_modify aevol_propagate



# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

aevol_run_SOURCES = aevol_run.cpp
aevol_run_CPPFLAGS = $(AM_CPPFLAGS)

aevol_run_IN2P3_SOURCES = aevol_run.cpp
aevol_run_IN2P3_CPPFLAGS = $(AM_CPPFLAGS)

aevol_run_debug_SOURCES = aevol_run.cpp
aevol_run_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_run_X11_SOURCES = aevol_run.cpp
aevol_run_X11_CPPFLAGS = $(AM_CPPFLAGS)

aevol_run_X11_debug_SOURCES = aevol_run.cpp
aevol_run_X11_debug_CPPFLAGS = $(AM_CPPFLAGS)


aevol_create_SOURCES = aevol_create.cpp
aevol_create_CPPFLAGS = $(AM_CPPFLAGS)


aevol_modify_SOURCES = aevol_modify.cpp
aevol_modify_CPPFLAGS = $(AM_CPPFLAGS)


aevol_propagate_SOURCES = aevol_propagate.cpp
aevol_propagate_CPPFLAGS = $(AM_CPPFLAGS)

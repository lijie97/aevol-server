############################################
#                 Variables                #
############################################
SUBDIRS = libaevol post_treatments

AM_CPPFLAGS = -I./
AM_CPPFLAGS += -I$(top_srcdir)/src/libaevol

############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS




############################################
# Set programs to build and library to use #
############################################


##### The below solution is no longer used since it generates a dependency problem #####
#AM_LDFLAGS = -L$(top_srcdir)/src/libaevol
#if WITH_RAEVOL
#  LDADD = -lraevol
#else # WITHOUT_RAEVOL
#  LDADD = -laevol
#endif # RAEVOL
#####  #####

if WITH_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libraevol.a
else # WITHOUT_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libaevol.a $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a
endif # RAEVOL


if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = raevol_X11_debug
else # WITHOUT_DEBUG
  bin_PROGRAMS = raevol_X11
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = raevol_debug
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = raevol_IN2P3
else # WITHOUT_IN2P3
  bin_PROGRAMS = raevol
endif # IN2P3
endif # DEBUG
endif # X
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = aevol_X11_debug
else # WITHOUT_DEBUG
  bin_PROGRAMS = aevol_X11
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = aevol_debug
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = aevol_IN2P3
else # WITHOUT_IN2P3
  bin_PROGRAMS = aevol
endif # IN2P3
endif # DEBUG
endif # X
endif # RAEVOL

bin_PROGRAMS += create_exp change_pop change_env change_exp_setup

CLEANFILES = aevol aevol_IN2P3 raevol raevol_IN2P3 aevol_debug raevol_debug aevol_X11 raevol_X11 aevol_X11_debug raevol_X11_debug
CLEANFILES += create_exp change_pop change_env change_exp_setup




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

aevoldir = $(bindir)
aevol_SOURCES = aevol.cpp
aevol_CPPFLAGS = $(AM_CPPFLAGS)

aevol_IN2P3dir = $(bindir)
aevol_IN2P3_SOURCES = aevol.cpp
aevol_IN2P3_CPPFLAGS = $(AM_CPPFLAGS)

aevol_debugdir = $(bindir)
aevol_debug_SOURCES = aevol.cpp
aevol_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_X11dir = $(bindir)
aevol_X11_SOURCES = aevol.cpp
aevol_X11_CPPFLAGS = $(AM_CPPFLAGS)

aevol_X11_debugdir = $(bindir)
aevol_X11_debug_SOURCES = aevol.cpp
aevol_X11_debug_CPPFLAGS = $(AM_CPPFLAGS)

raevoldir = $(bindir)
raevol_SOURCES = aevol.cpp
raevol_CPPFLAGS = $(AM_CPPFLAGS)

raevol_IN2P3dir = $(bindir)
raevol_IN2P3_SOURCES = aevol.cpp
raevol_IN2P3_CPPFLAGS = $(AM_CPPFLAGS)

raevol_debugdir = $(bindir)
raevol_debug_SOURCES = aevol.cpp
raevol_debug_CPPFLAGS = $(AM_CPPFLAGS)

raevol_X11dir = $(bindir)
raevol_X11_SOURCES = aevol.cpp
raevol_X11_CPPFLAGS = $(AM_CPPFLAGS)

raevol_X11_debugdir = $(bindir)
raevol_X11_debug_SOURCES = aevol.cpp
raevol_X11_debug_CPPFLAGS = $(AM_CPPFLAGS)


create_exp_dir = $(bindir)
create_exp_SOURCES = create_exp.cpp
create_exp_CPPFLAGS = $(AM_CPPFLAGS)

change_pop_dir = $(bindir)
change_pop_SOURCES = change_pop.cpp
change_pop_CPPFLAGS = $(AM_CPPFLAGS)

change_env_dir = $(bindir)
change_env_SOURCES = change_env.cpp
change_env_CPPFLAGS = $(AM_CPPFLAGS)

change_exp_setup_dir = $(bindir)
change_exp_setup_SOURCES = change_exp_setup.cpp
change_exp_setup_CPPFLAGS = $(AM_CPPFLAGS)

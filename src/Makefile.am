######################    variables    ########################
SUBDIRS = lib-aevol post_treatments

AM_CPPFLAGS = -I./
AM_CPPFLAGS += -I$(top_srcdir)/src/lib-aevol

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
if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = raevol_X11_debug
  LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = raevol_X11
  LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11.a
endif # DEBUG
else # WITHOUT_X11
if WITH_DEBUG
  bin_PROGRAMS = raevol_debug
  LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = raevol_IN2P3
  LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = raevol
  LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol.a
endif # IN2P3
endif # DEBUG
endif # X
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = aevol_X11_debug
  LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = aevol_X11
  LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11.a
endif # DEBUG
else # WITHOUT_X11
if WITH_DEBUG
  bin_PROGRAMS = aevol_debug
  LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = aevol_IN2P3
  LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = aevol
  LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol.a
endif # IN2P3
endif # DEBUG
endif # X
endif # RAEVOL

CLEANFILES = aevol aevol_IN2P3 raevol raevol_IN2P3 aevol_debug raevol_debug aevol_X11 raevol_X11 aevol_X11_debug raevol_X11_debug





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

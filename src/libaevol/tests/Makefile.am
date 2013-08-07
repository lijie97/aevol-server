if WITH_CPPUNIT

############################################
#         Set aevol library to use         #
############################################


AM_CPPFLAGS = -I./
AM_CPPFLAGS += -I$(top_srcdir)/src/libaevol $(CPPUNIT_CFLAGS)


AM_LDFLAGS = -L$(top_srcdir)/src/libaevol 

if WITH_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libraevol.a
else # WITHOUT_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libaevol.a
endif # RAEVOL

LDADD += $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a $(CPPUNIT_LIBS)


############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS


############################################
#           Set programs to build          #
############################################
bin_PROGRAMS = run_tests
CLEANFILES = $(bin_PROGRAMS)
TESTS = $(bin_PROGRAMS)



# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

run_testsdir = $(bindir)
run_tests_SOURCES = run_tests.cpp Test_ae_jumping_mt.cpp Test_ae_jumping_mt.h
run_tests_CPPFLAGS = $(AM_CPPFLAGS)


endif

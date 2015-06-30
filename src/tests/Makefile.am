############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)
gtest_dir = gtest-1.7.0

SUBDIRS = $(gtest_dir)


# We define the global AM_CPPFLAGS and not target-specific flags because
# everything we compile includes from these directories.
AM_CPPFLAGS += -I$(top_srcdir)/src/libaevol
AM_CPPFLAGS += -isystem $(gtest_dir) -isystem $(gtest_dir)/include

# Modifies compiler and linker flags for pthreads compatibility.
if HAVE_PTHREADS
  AM_CXXFLAGS += @PTHREAD_CFLAGS@ -DGTEST_HAS_PTHREAD=1
  AM_LIBS = @PTHREAD_LIBS@
else
  AM_CXXFLAGS += -DGTEST_HAS_PTHREAD=0
endif

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
#           Set programs to build          #
############################################
check_PROGRAMS = sampleTest IndividualTest TranscriptionTranslationTest
CLEANFILES = $(check_PROGRAMS)
TESTS = $(check_PROGRAMS)


sampleTest_SOURCES = sampleTest.cpp
sampleTest_LDADD = $(gtest_dir)/lib/libgtest_main.a \
                   $(gtest_dir)/lib/libgtest.a

IndividualTest_SOURCES = IndividualTest.cpp
IndividualTest_LDADD =  $(gtest_dir)/lib/libgtest_main.a \
                        $(gtest_dir)/lib/libgtest.a \
                        $(top_srcdir)/src/libaevol/libaevol.a \
                        $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a

TranscriptionTranslationTest_SOURCES = TranscriptionTranslationTest.cpp
TranscriptionTranslationTest_LDADD = $(gtest_dir)/lib/libgtest_main.a \
			                         $(gtest_dir)/lib/libgtest.a \
			                         $(top_srcdir)/src/libaevol/libaevol.a \
			                         $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a

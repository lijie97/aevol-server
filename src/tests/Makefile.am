############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)

SUBDIRS = gtest-1.7.0


# We define the global AM_CPPFLAGS and not target-specific flags because
# everything we compile includes from these directories.
AM_CPPFLAGS += -I$(srcdir) -I$(srcdir)/include

# Modifies compiler and linker flags for pthreads compatibility.
if HAVE_PTHREADS
  AM_CXXFLAGS += @PTHREAD_CFLAGS@ -DGTEST_HAS_PTHREAD=1
  AM_LIBS = @PTHREAD_LIBS@
else
  AM_CXXFLAGS += -DGTEST_HAS_PTHREAD=0
endif




############################################
#           Set programs to build          #
############################################
check_PROGRAMS = sampleTest
CLEANFILES = $(check_PROGRAMS)
TESTS = $(check_PROGRAMS)


sampleTest_SOURCES = sampleTest.cpp
sampleTest_LDADD = gtest-1.7.0/lib/libgtest_main.a \
                   gtest-1.7.0/lib/libgtest.a

############################################
#         Set aevol library to use         #
############################################
AM_LDFLAGS = -L$(top_srcdir)/src/libaevol

AM_CPPFLAGS = -I./
AM_CPPFLAGS += -I$(top_srcdir)/src/libaevol

if WITH_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libraevol.a
else # WITHOUT_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libaevol.a $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a
endif # RAEVOL


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
if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug rextract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug rextract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol.a
endif # IN2P3
endif # DEBUG
endif # X
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug extract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug extract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol.a
endif # IN2P3
endif # DEBUG
endif # X
endif # RAEVOL

CLEANFILES = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug
CLEANFILES += rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats
CLEANFILES += robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug extract_debug
CLEANFILES += robustness ancstats create_eps fixed_mutations lineage compute_pop_stats extract




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

ancstats_SOURCES = ancstats.cpp
ancstats_CPPFLAGS = $(AM_CPPFLAGS)

ancstats_debug_SOURCES = ancstats.cpp
ancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)

rancstats_SOURCES = ancstats.cpp
rancstats_CPPFLAGS = $(AM_CPPFLAGS)

rancstats_debug_SOURCES = ancstats.cpp
rancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)


create_eps_SOURCES = create_eps.cpp
create_eps_CPPFLAGS = $(AM_CPPFLAGS)

create_eps_debug_SOURCES = create_eps.cpp
create_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)

rcreate_eps_SOURCES = create_eps.cpp
rcreate_eps_CPPFLAGS = $(AM_CPPFLAGS)

rcreate_eps_debug_SOURCES = create_eps.cpp
rcreate_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)


extract_SOURCES = extract.cpp
extract_CPPFLAGS = $(AM_CPPFLAGS)

extract_debug_SOURCES = extract.cpp
extract_debug_CPPFLAGS = $(AM_CPPFLAGS)

rextract_SOURCES = extract.cpp
rextract_CPPFLAGS = $(AM_CPPFLAGS)

rextract_debug_SOURCES = extract.cpp
rextract_debug_CPPFLAGS = $(AM_CPPFLAGS)


fixed_mutations_SOURCES = fixed_mutations.cpp
fixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

fixed_mutations_debug_SOURCES = fixed_mutations.cpp
fixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)

rfixed_mutations_SOURCES = fixed_mutations.cpp
rfixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

rfixed_mutations_debug_SOURCES = fixed_mutations.cpp
rfixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)

lineage_SOURCES = lineage.cpp
lineage_CPPFLAGS = $(AM_CPPFLAGS)

lineage_debug_SOURCES = lineage.cpp
lineage_debug_CPPFLAGS = $(AM_CPPFLAGS)

rlineage_SOURCES = lineage.cpp
rlineage_CPPFLAGS = $(AM_CPPFLAGS)

rlineage_debug_SOURCES = lineage.cpp
rlineage_debug_CPPFLAGS = $(AM_CPPFLAGS)


#compute_pop_stats_HEADERS = population_statistics.h
compute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
compute_pop_stats_SOURCES += population_statistics.h
compute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

#compute_pop_stats_debug_HEADERS = population_statistics.h
compute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
compute_pop_stats_debug_SOURCES += population_statistics.h
compute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rcompute_pop_stats_HEADERS = population_statistics.h
rcompute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
rcompute_pop_stats_SOURCES += population_statistics.h
rcompute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

#rcompute_pop_stats_debug_HEADERS = population_statistics.h
rcompute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
rcompute_pop_stats_debug_SOURCES += population_statistics.h
rcompute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)


#robustness_HEADERS = population_statistics.h
robustness_SOURCES = population_statistics.cpp robustness.cpp
robustness_SOURCES +=population_statistics.h
robustness_CPPFLAGS = $(AM_CPPFLAGS)

#robustness_debug_HEADERS = population_statistics.h
robustness_debug_SOURCES = population_statistics.cpp robustness.cpp
robustness_debug_SOURCES += population_statistics.h
robustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rrobustness_HEADERS = population_statistics.h
rrobustness_SOURCES = population_statistics.cpp robustness.cpp
rrobustness_SOURCES += population_statistics.h
rrobustness_CPPFLAGS = $(AM_CPPFLAGS)

#rrobustness_debug_HEADERS = population_statistics.h
rrobustness_debug_SOURCES = population_statistics.cpp robustness.cpp
rrobustness_debug_SOURCES += population_statistics.h
rrobustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

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


####################################
# Add hardening option to compiler #
####################################
# This replaces insecure unlimited length buffer function calls with
# length-limited ones (required for debian and not a bad idea anyway).
AM_CPPFLAGS += -D_FORTIFY_SOURCE=2


############################################
#           Set programs to build          #
############################################
if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  libexec_PROGRAMS = aevol_misc_rrobustness_debug aevol_misc_rancstats_debug aevol_misc_rcreate_eps_debug aevol_misc_rfixed_mutations_debug aevol_misc_rlineage_debug aevol_misc_rcompute_pop_stats_debug aevol_misc_rextract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11-debug.a
else # WITHOUT_DEBUG
  libexec_PROGRAMS = aevol_misc_rrobustness aevol_misc_rancstats aevol_misc_rcreate_eps aevol_misc_rfixed_mutations aevol_misc_rlineage aevol_misc_rcompute_pop_stats aevol_misc_rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  libexec_PROGRAMS = aevol_misc_rrobustness_debug aevol_misc_rancstats_debug aevol_misc_rcreate_eps_debug aevol_misc_rfixed_mutations_debug aevol_misc_rlineage_debug aevol_misc_rcompute_pop_stats_debug aevol_misc_rextract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  libexec_PROGRAMS = aevol_misc_rrobustness aevol_misc_rancstats aevol_misc_rcreate_eps aevol_misc_rfixed_mutations aevol_misc_rlineage aevol_misc_rcompute_pop_stats aevol_misc_rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-IN2P3.a
else # WITHOUT_IN2P3
  libexec_PROGRAMS = aevol_misc_rrobustness aevol_misc_rancstats aevol_misc_rcreate_eps aevol_misc_rfixed_mutations aevol_misc_rlineage aevol_misc_rcompute_pop_stats aevol_misc_rextract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol.a
endif # IN2P3
endif # DEBUG
endif # X
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  libexec_PROGRAMS = aevol_misc_robustness_debug aevol_misc_ancstats_debug aevol_misc_create_eps_debug aevol_misc_fixed_mutations_debug aevol_misc_lineage_debug aevol_misc_compute_pop_stats_debug aevol_misc_extract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11-debug.a
else # WITHOUT_DEBUG
  libexec_PROGRAMS = aevol_misc_robustness aevol_misc_ancstats aevol_misc_create_eps aevol_misc_fixed_mutations aevol_misc_lineage aevol_misc_compute_pop_stats aevol_misc_extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  libexec_PROGRAMS = aevol_misc_robustness_debug aevol_misc_ancstats_debug aevol_misc_create_eps_debug aevol_misc_fixed_mutations_debug aevol_misc_lineage_debug aevol_misc_compute_pop_stats_debug aevol_misc_extract_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  libexec_PROGRAMS = aevol_misc_robustness aevol_misc_ancstats aevol_misc_create_eps aevol_misc_fixed_mutations aevol_misc_lineage aevol_misc_compute_pop_stats aevol_misc_extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-IN2P3.a
else # WITHOUT_IN2P3
  libexec_PROGRAMS = aevol_misc_robustness aevol_misc_ancstats aevol_misc_create_eps aevol_misc_fixed_mutations aevol_misc_lineage aevol_misc_compute_pop_stats aevol_misc_extract
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol.a
endif # IN2P3
endif # DEBUG
endif # X
endif # RAEVOL

CLEANFILES = aevol_misc_rrobustness_debug aevol_misc_rancstats_debug aevol_misc_rcreate_eps_debug aevol_misc_rfixed_mutations_debug aevol_misc_rlineage_debug aevol_misc_rcompute_pop_stats_debug
CLEANFILES += aevol_misc_rrobustness aevol_misc_rancstats aevol_misc_rcreate_eps aevol_misc_rfixed_mutations aevol_misc_rlineage aevol_misc_rcompute_pop_stats
CLEANFILES += aevol_misc_robustness_debug aevol_misc_ancstats_debug aevol_misc_create_eps_debug aevol_misc_fixed_mutations_debug aevol_misc_lineage_debug compute_pop_stats_debug extract_debug
CLEANFILES += aevol_misc_robustness aevol_misc_ancstats aevol_misc_create_eps aevol_misc_fixed_mutations aevol_misc_lineage aevol_misc_compute_pop_stats aevol_misc_extract




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

aevol_misc_ancstats_SOURCES = ancstats.cpp
aevol_misc_ancstats_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_ancstats_debug_SOURCES = ancstats.cpp
aevol_misc_ancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rancstats_SOURCES = ancstats.cpp
aevol_misc_rancstats_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rancstats_debug_SOURCES = ancstats.cpp
aevol_misc_rancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)


aevol_misc_create_eps_SOURCES = create_eps.cpp
aevol_misc_create_eps_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_create_eps_debug_SOURCES = create_eps.cpp
aevol_misc_create_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rcreate_eps_SOURCES = create_eps.cpp
aevol_misc_rcreate_eps_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rcreate_eps_debug_SOURCES = create_eps.cpp
aevol_misc_rcreate_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)


aevol_misc_extract_SOURCES = extract.cpp
aevol_misc_extract_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_extract_debug_SOURCES = extract.cpp
aevol_misc_extract_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rextract_SOURCES = extract.cpp
aevol_misc_rextract_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rextract_debug_SOURCES = extract.cpp
aevol_misc_rextract_debug_CPPFLAGS = $(AM_CPPFLAGS)


aevol_misc_fixed_mutations_SOURCES = fixed_mutations.cpp
aevol_misc_fixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_fixed_mutations_debug_SOURCES = fixed_mutations.cpp
aevol_misc_fixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rfixed_mutations_SOURCES = fixed_mutations.cpp
aevol_misc_rfixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rfixed_mutations_debug_SOURCES = fixed_mutations.cpp
aevol_misc_rfixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_lineage_SOURCES = lineage.cpp
aevol_misc_lineage_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_lineage_debug_SOURCES = lineage.cpp
aevol_misc_lineage_debug_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rlineage_SOURCES = lineage.cpp
aevol_misc_rlineage_CPPFLAGS = $(AM_CPPFLAGS)

aevol_misc_rlineage_debug_SOURCES = lineage.cpp
aevol_misc_rlineage_debug_CPPFLAGS = $(AM_CPPFLAGS)


#aevol_misc_compute_pop_stats_HEADERS = population_statistics.h
aevol_misc_compute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
aevol_misc_compute_pop_stats_SOURCES += population_statistics.h
aevol_misc_compute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_compute_pop_stats_debug_HEADERS = population_statistics.h
aevol_misc_compute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
aevol_misc_compute_pop_stats_debug_SOURCES += population_statistics.h
aevol_misc_compute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_rcompute_pop_stats_HEADERS = population_statistics.h
aevol_misc_rcompute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
aevol_misc_rcompute_pop_stats_SOURCES += population_statistics.h
aevol_misc_rcompute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_rcompute_pop_stats_debug_HEADERS = population_statistics.h
aevol_misc_rcompute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
aevol_misc_rcompute_pop_stats_debug_SOURCES += population_statistics.h
aevol_misc_rcompute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)


#aevol_misc_robustness_HEADERS = population_statistics.h
aevol_misc_robustness_SOURCES = population_statistics.cpp robustness.cpp
aevol_misc_robustness_SOURCES +=population_statistics.h
aevol_misc_robustness_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_robustness_debug_HEADERS = population_statistics.h
aevol_misc_robustness_debug_SOURCES = population_statistics.cpp robustness.cpp
aevol_misc_robustness_debug_SOURCES += population_statistics.h
aevol_misc_robustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_rrobustness_HEADERS = population_statistics.h
aevol_misc_rrobustness_SOURCES = population_statistics.cpp robustness.cpp
aevol_misc_rrobustness_SOURCES += population_statistics.h
aevol_misc_rrobustness_CPPFLAGS = $(AM_CPPFLAGS)

#aevol_misc_rrobustness_debug_HEADERS = population_statistics.h
aevol_misc_rrobustness_debug_SOURCES = population_statistics.cpp robustness.cpp
aevol_misc_rrobustness_debug_SOURCES += population_statistics.h
aevol_misc_rrobustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

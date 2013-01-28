############################################
# Set programs to build and library to use #
############################################
AM_LDFLAGS = -L$(top_srcdir)/src/libaevol

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

if WITH_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libraevol.a
else # WITHOUT_RAEVOL
  LDADD = $(top_srcdir)/src/libaevol/libaevol.a $(top_srcdir)/src/libaevol/SFMT-src-1.4/libsfmt.a
endif # RAEVOL

#if WITH_RAEVOL
#if WITH_X
#if WITH_DEBUG
#  bin_PROGRAMS =  rancstats_debug rcreate_eps_debug rextract_debug rfixed_mutations_debug rgenome_stats_debug rget_best_indiv_debug
#  bin_PROGRAMS += rlineage_debug rmutagenesis_debug rrobustness_debug rset_genome_from_indiv_debug rview_generation_debug rcompute_pop_stats_debug
#else # WITHOUT_DEBUG
#  bin_PROGRAMS = rancstats rcreate_eps rextract rfixed_mutations rgenome_stats rget_best_indiv rlineage rmutagenesis rrobustness rset_genome_from_indiv rview_generation rcompute_pop_stats
#endif # DEBUG
#else # WITHOUT_X
#if WITH_DEBUG
#  bin_PROGRAMS =  rancstats_debug rcreate_eps_debug rextract_debug rfixed_mutations_debug rgenome_stats_debug rget_best_indiv_debug
#  bin_PROGRAMS += rlineage_debug rmutagenesis_debug rrobustness_debug rset_genome_from_indiv_debug rcompute_pop_stats_debug
#else # WITHOUT_DEBUG
#if WITH_IN2P3
#  bin_PROGRAMS = rancstats rcreate_eps rextract rfixed_mutations rgenome_stats rget_best_indiv rlineage rmutagenesis rrobustness rset_genome_from_indiv rcompute_pop_stats
#else # WITHOUT_IN2P3
#  bin_PROGRAMS = rancstats rcreate_eps rextract rfixed_mutations rgenome_stats rget_best_indiv rlineage rmutagenesis rrobustness rset_genome_from_indiv rcompute_pop_stats
#endif # IN2P3
#endif # DEBUG
#endif # X
#else # WITHOUT_RAEVOL
#if WITH_X
#if WITH_DEBUG
#  bin_PROGRAMS =  ancstats_debug create_eps_debug extract_debug fixed_mutations_debug genome_stats_debug get_best_indiv_debug
#  bin_PROGRAMS += lineage_debug mutagenesis_debug mutationalrobustness_debug robustness_debug set_genome_from_indiv_debug view_generation_debug compute_pop_stats_debug
#else # WITHOUT_DEBUG
#  bin_PROGRAMS = ancstats create_eps extract fixed_mutations genome_stats get_best_indiv lineage mutagenesis mutationalrobustness robustness set_genome_from_indiv view_generation compute_pop_stats
#endif # DEBUG
#else # WITHOUT_X
#if WITH_DEBUG
#  bin_PROGRAMS =  ancstats_debug create_eps_debug extract_debug fixed_mutations_debug genome_stats_debug get_best_indiv_debug
#  bin_PROGRAMS += lineage_debug mutagenesis_debug mutationalrobustness_debug robustness_debug set_genome_from_indiv_debug compute_pop_stats_debug
#else # WITHOUT_DEBUG
#if WITH_IN2P3
#  bin_PROGRAMS = ancstats create_eps extract fixed_mutations genome_stats get_best_indiv lineage mutagenesis mutationalrobustness robustness set_genome_from_indiv compute_pop_stats
#else # WITHOUT_IN2P3
#  bin_PROGRAMS = ancstats create_eps extract fixed_mutations genome_stats get_best_indiv lineage mutagenesis mutationalrobustness robustness set_genome_from_indiv compute_pop_stats
#endif # IN2P3
#endif # DEBUG
#endif # X
#endif # RAEVOL

if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-raevol.a
endif # IN2P3
endif # DEBUG
endif # X
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  bin_PROGRAMS = robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11-debug.a
else # WITHOUT_DEBUG
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-X11.a
endif # DEBUG
else # WITHOUT_X
if WITH_DEBUG
  bin_PROGRAMS = robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol-IN2P3.a
else # WITHOUT_IN2P3
  bin_PROGRAMS = robustness ancstats create_eps fixed_mutations lineage compute_pop_stats
#	LDADD = $(top_srcdir)/src/lib-aevol/lib-aevol.a
endif # IN2P3
endif # DEBUG
endif # X
endif # RAEVOL

CLEANFILES = rrobustness_debug rancstats_debug rcreate_eps_debug rfixed_mutations_debug rlineage_debug rcompute_pop_stats_debug
CLEANFILES += rrobustness rancstats rcreate_eps rfixed_mutations rlineage rcompute_pop_stats
CLEANFILES += robustness_debug ancstats_debug create_eps_debug fixed_mutations_debug lineage_debug compute_pop_stats_debug
CLEANFILES += robustness ancstats create_eps fixed_mutations lineage compute_pop_stats




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

ancstatsdir = $(bindir)
ancstats_SOURCES = ancstats.cpp
ancstats_CPPFLAGS = $(AM_CPPFLAGS)

ancstats_debugdir = $(bindir)
ancstats_debug_SOURCES = ancstats.cpp
ancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)

rancstatsdir = $(bindir)
rancstats_SOURCES = ancstats.cpp
rancstats_CPPFLAGS = $(AM_CPPFLAGS)

rancstats_debugdir = $(bindir)
rancstats_debug_SOURCES = ancstats.cpp
rancstats_debug_CPPFLAGS = $(AM_CPPFLAGS)


create_epsdir = $(bindir)
create_eps_SOURCES = create_eps.cpp
create_eps_CPPFLAGS = $(AM_CPPFLAGS)

create_eps_debugdir = $(bindir)
create_eps_debug_SOURCES = create_eps.cpp
create_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)

rcreate_epsdir = $(bindir)
rcreate_eps_SOURCES = create_eps.cpp
rcreate_eps_CPPFLAGS = $(AM_CPPFLAGS)

rcreate_eps_debugdir = $(bindir)
rcreate_eps_debug_SOURCES = create_eps.cpp
rcreate_eps_debug_CPPFLAGS = $(AM_CPPFLAGS)


#extractdir = $(bindir)
#extract_SOURCES = extract.cpp
#extract_CPPFLAGS = $(AM_CPPFLAGS)

#extract_debugdir = $(bindir)
#extract_debug_SOURCES = extract.cpp
#extract_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rextractdir = $(bindir)
#rextract_SOURCES = extract.cpp
#rextract_CPPFLAGS = $(AM_CPPFLAGS)

#rextract_debugdir = $(bindir)
#rextract_debug_SOURCES = extract.cpp
#rextract_debug_CPPFLAGS = $(AM_CPPFLAGS)


fixed_mutationsdir = $(bindir)
fixed_mutations_SOURCES = fixed_mutations.cpp
fixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

fixed_mutations_debugdir = $(bindir)
fixed_mutations_debug_SOURCES = fixed_mutations.cpp
fixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)

rfixed_mutationsdir = $(bindir)
rfixed_mutations_SOURCES = fixed_mutations.cpp
rfixed_mutations_CPPFLAGS = $(AM_CPPFLAGS)

rfixed_mutations_debugdir = $(bindir)
rfixed_mutations_debug_SOURCES = fixed_mutations.cpp
rfixed_mutations_debug_CPPFLAGS = $(AM_CPPFLAGS)


#genome_statsdir = $(bindir)
#genome_stats_SOURCES = genome_stats.cpp
#genome_stats_CPPFLAGS = $(AM_CPPFLAGS)

#genome_stats_debugdir = $(bindir)
#genome_stats_debug_SOURCES = genome_stats.cpp
#genome_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rgenome_statsdir = $(bindir)
#rgenome_stats_SOURCES = genome_stats.cpp
#rgenome_stats_CPPFLAGS = $(AM_CPPFLAGS)

#rgenome_stats_debugdir = $(bindir)
#rgenome_stats_debug_SOURCES = genome_stats.cpp
#rgenome_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)



#get_best_indivdir = $(bindir)
#get_best_indiv_SOURCES = get_best_indiv.cpp
#get_best_indiv_CPPFLAGS = $(AM_CPPFLAGS)

#get_best_indiv_debugdir = $(bindir)
#get_best_indiv_debug_SOURCES = get_best_indiv.cpp
#get_best_indiv_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rget_best_indivdir = $(bindir)
#rget_best_indiv_SOURCES = get_best_indiv.cpp
#rget_best_indiv_CPPFLAGS = $(AM_CPPFLAGS)

#rget_best_indiv_debugdir = $(bindir)
#rget_best_indiv_debug_SOURCES = get_best_indiv.cpp
#rget_best_indiv_debug_CPPFLAGS = $(AM_CPPFLAGS)


lineagedir = $(bindir)
lineage_SOURCES = lineage.cpp
lineage_CPPFLAGS = $(AM_CPPFLAGS)

lineage_debugdir = $(bindir)
lineage_debug_SOURCES = lineage.cpp
lineage_debug_CPPFLAGS = $(AM_CPPFLAGS)

rlineagedir = $(bindir)
rlineage_SOURCES = lineage.cpp
rlineage_CPPFLAGS = $(AM_CPPFLAGS)

rlineage_debugdir = $(bindir)
rlineage_debug_SOURCES = lineage.cpp
rlineage_debug_CPPFLAGS = $(AM_CPPFLAGS)


#mutagenesisdir = $(bindir)
#mutagenesis_SOURCES = mutagenesis.cpp
#mutagenesis_CPPFLAGS = $(AM_CPPFLAGS)

#mutagenesis_debugdir = $(bindir)
#mutagenesis_debug_SOURCES = mutagenesis.cpp
#mutagenesis_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rmutagenesisdir = $(bindir)
#rmutagenesis_SOURCES = mutagenesis.cpp
#rmutagenesis_CPPFLAGS = $(AM_CPPFLAGS)

#rmutagenesis_debugdir = $(bindir)
#rmutagenesis_debug_SOURCES = mutagenesis.cpp
#rmutagenesis_debug_CPPFLAGS = $(AM_CPPFLAGS)


#mutationalrobustnessdir = $(bindir)
#mutationalrobustness_SOURCES = mutationalrobustness.cpp
#mutationalrobustness_CPPFLAGS = $(AM_CPPFLAGS)

#mutationalrobustness_debugdir = $(bindir)
#mutationalrobustness_debug_SOURCES = mutationalrobustness.cpp
#mutationalrobustness_debug_CPPFLAGS = $(AM_CPPFLAGS)


#robustnessdir = $(bindir)
#robustness_SOURCES = robustness.cpp
#robustness_CPPFLAGS = $(AM_CPPFLAGS)

#robustness_debugdir = $(bindir)
#robustness_debug_SOURCES = robustness.cpp
#robustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rrobustnessdir = $(bindir)
#rrobustness_SOURCES = robustness.cpp
#rrobustness_CPPFLAGS = $(AM_CPPFLAGS)

#rrobustness_debugdir = $(bindir)
#rrobustness_debug_SOURCES = robustness.cpp
#rrobustness_debug_CPPFLAGS = $(AM_CPPFLAGS)


#set_genome_from_indivdir = $(bindir)
#set_genome_from_indiv_SOURCES = set_genome_from_indiv.cpp
#set_genome_from_indiv_CPPFLAGS = $(AM_CPPFLAGS)

#set_genome_from_indiv_debugdir = $(bindir)
#set_genome_from_indiv_debug_SOURCES = set_genome_from_indiv.cpp
#set_genome_from_indiv_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rset_genome_from_indivdir = $(bindir)
#rset_genome_from_indiv_SOURCES = set_genome_from_indiv.cpp
#rset_genome_from_indiv_CPPFLAGS = $(AM_CPPFLAGS)

#rset_genome_from_indiv_debugdir = $(bindir)
#rset_genome_from_indiv_debug_SOURCES = set_genome_from_indiv.cpp
#rset_genome_from_indiv_debug_CPPFLAGS = $(AM_CPPFLAGS)


#view_generationdir = $(bindir)
#view_generation_SOURCES = view_generation.cpp
#view_generation_CPPFLAGS = $(AM_CPPFLAGS)

#view_generation_debugdir = $(bindir)
#view_generation_debug_SOURCES = view_generation.cpp
#view_generation_debug_CPPFLAGS = $(AM_CPPFLAGS)

#rview_generationdir = $(bindir)
#rview_generation_SOURCES = view_generation.cpp
#rview_generation_CPPFLAGS = $(AM_CPPFLAGS)

#rview_generation_debugdir = $(bindir)
#rview_generation_debug_SOURCES = view_generation.cpp
#rview_generation_debug_CPPFLAGS = $(AM_CPPFLAGS)


compute_pop_statsdir = $(bindir)
compute_pop_stats_HEADERS = population_statistics.h
compute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
compute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

compute_pop_stats_debugdir = $(bindir)
compute_pop_stats_debug_HEADERS = population_statistics.h
compute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
compute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)

rcompute_pop_statsdir = $(bindir)
rcompute_pop_stats_HEADERS = population_statistics.h
rcompute_pop_stats_SOURCES = population_statistics.cpp compute_pop_stats.cpp
rcompute_pop_stats_CPPFLAGS = $(AM_CPPFLAGS)

rcompute_pop_stats_debugdir = $(bindir)
rcompute_pop_stats_debug_HEADERS = population_statistics.h
rcompute_pop_stats_debug_SOURCES = population_statistics.cpp compute_pop_stats.cpp
rcompute_pop_stats_debug_CPPFLAGS = $(AM_CPPFLAGS)


robustnessdir = $(bindir)
robustness_HEADERS = population_statistics.h
robustness_SOURCES = population_statistics.cpp robustness.cpp
robustness_CPPFLAGS = $(AM_CPPFLAGS)

robustness_debugdir = $(bindir)
robustness_debug_HEADERS = population_statistics.h
robustness_debug_SOURCES = population_statistics.cpp robustness.cpp
robustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

rrobustnessdir = $(bindir)
rrobustness_HEADERS = population_statistics.h
rrobustness_SOURCES = population_statistics.cpp robustness.cpp
rrobustness_CPPFLAGS = $(AM_CPPFLAGS)

rrobustness_debugdir = $(bindir)
rrobustness_debug_HEADERS = population_statistics.h
rrobustness_debug_SOURCES = population_statistics.cpp robustness.cpp
rrobustness_debug_CPPFLAGS = $(AM_CPPFLAGS)

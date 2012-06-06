########################
# Set library to build #
########################
if WITH_RAEVOL
if WITH_X
if WITH_DEBUG
  noinst_LIBRARIES = lib-raevol-X11-debug.a
else # WITHOUT_DEBUG
  noinst_LIBRARIES = lib-raevol-X11.a
endif # DEBUG
else # WITHOUT_X11
if WITH_DEBUG
  noinst_LIBRARIES = lib-raevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  noinst_LIBRARIES = lib-raevol-IN2P3.a
else # WITHOUT_IN2P3
  noinst_LIBRARIES = lib-raevol.a
endif # IN2P3
endif # DEBUG
endif # X11
else # WITHOUT_RAEVOL
if WITH_X
if WITH_DEBUG
  noinst_LIBRARIES = lib-aevol-X11-debug.a
else # WITHOUT_DEBUG
  noinst_LIBRARIES = lib-aevol-X11.a
endif # DEBUG
else # WITHOUT_X11
if WITH_DEBUG
  noinst_LIBRARIES = lib-aevol-debug.a
else # WITHOUT_DEBUG
if WITH_IN2P3
  noinst_LIBRARIES = lib-aevol-IN2P3.a
else # WITHOUT_IN2P3
  noinst_LIBRARIES = lib-aevol.a
endif # IN2P3
endif # DEBUG
endif # X11
endif # RAEVOL


AM_CPPFLAGS = -I./

############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS



common_headers = ae_align.h ae_codon.h ae_common.h ae_dna.h ae_dump.h ae_enums.h ae_environment.h ae_env_segment.h ae_fuzzy_set.h ae_gaussian.h ae_genetic_unit.h ae_grid_cell.h ae_individual.h ae_list.h ae_logs.h ae_macros.h ae_mutation.h ae_object.h ae_params.h ae_params_init.h ae_params_record.h ae_param_loader.h ae_param_overloader.h ae_phenotype.h ae_population.h ae_point_1d.h ae_point_2d.h ae_protein.h ae_rand_mt.h ae_replication_report.h ae_dna_replic_report.h ae_rna.h ae_simulation.h ae_stats.h ae_stat_record.h ae_string.h ae_tree.h ae_utils.h ae_vis_a_vis.h

common_sources = ae_align.cpp ae_codon.cpp ae_common.cpp ae_dna.cpp ae_dump.cpp ae_environment.cpp ae_fuzzy_set.cpp ae_genetic_unit.cpp ae_grid_cell.cpp ae_individual.cpp ae_logs.cpp ae_mutation.cpp ae_params.cpp ae_params_init.cpp ae_params_record.cpp ae_param_loader.cpp ae_param_overloader.cpp ae_phenotype.cpp ae_population.cpp ae_protein.cpp ae_rand_mt.cpp ae_replication_report.cpp ae_dna_replic_report.cpp ae_rna.cpp ae_simulation.cpp ae_stats.cpp ae_stat_record.cpp ae_string.cpp ae_tree.cpp  ae_vis_a_vis.cpp


X11_headers = ae_fuzzy_set_X11.h ae_individual_X11.h ae_population_X11.h ae_simulation_X11.h ae_X11_window.h
X11_sources = ae_fuzzy_set_X11.cpp ae_individual_X11.cpp ae_population_X11.cpp ae_simulation_X11.cpp ae_X11_window.cpp

raevol_headers = ae_array_short.h ae_individual_R.h ae_influence_R.h ae_matrix_double.h ae_protein_R.h ae_rna_R.h
raevol_sources = ae_individual_R.cpp ae_influence_R.cpp ae_protein_R.cpp ae_rna_R.cpp


X11_regulation_headers = $(X11_headers) $(raevol_headers) ae_individual_R_X11.h
X11_regulation_sources = $(X11_sources) $(raevol_sources) ae_individual_R_X11.cpp




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

lib_aevol_adir = ./
lib_aevol_a_HEADERS = $(common_headers)
lib_aevol_a_SOURCES = $(common_sources)
lib_aevol_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_aevol_IN2P3_adir = ./
lib_aevol_IN2P3_a_HEADERS = $(common_headers)
lib_aevol_IN2P3_a_SOURCES = $(common_sources)
lib_aevol_IN2P3_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_aevol_debug_adir = ./
lib_aevol_debug_a_HEADERS = $(common_headers)
lib_aevol_debug_a_SOURCES = $(common_sources)
lib_aevol_debug_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_aevol_X11_adir = ./
lib_aevol_X11_a_HEADERS = $(common_headers) $(X11_headers)
lib_aevol_X11_a_SOURCES = $(common_sources) $(X11_sources)
lib_aevol_X11_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_aevol_X11_debug_adir = ./
lib_aevol_X11_debug_a_HEADERS = $(common_headers) $(X11_headers)
lib_aevol_X11_debug_a_SOURCES = $(common_sources) $(X11_sources)
lib_aevol_X11_debug_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_raevol_adir = ./
lib_raevol_a_HEADERS = $(common_headers) $(raevol_headers)
lib_raevol_a_SOURCES = $(common_sources) $(raevol_sources)
lib_raevol_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_raevol_IN2P3_adir = ./
lib_raevol_IN2P3_a_HEADERS = $(common_headers) $(raevol_headers)
lib_raevol_IN2P3_a_SOURCES = $(common_sources) $(raevol_sources)
lib_raevol_IN2P3_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_raevol_debug_adir = ./
lib_raevol_debug_a_HEADERS = $(common_headers) $(raevol_headers)
lib_raevol_debug_a_SOURCES = $(common_sources) $(raevol_sources)
lib_raevol_debug_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_raevol_X11_adir = ./
lib_raevol_X11_a_HEADERS = $(common_headers) $(X11_regulation_headers)
lib_raevol_X11_a_SOURCES = $(common_sources) $(X11_regulation_sources)
lib_raevol_X11_a_CPPFLAGS = $(AM_CPPFLAGS)

lib_raevol_X11_debug_adir = ./
lib_raevol_X11_debug_a_HEADERS = $(common_headers) $(X11_regulation_headers)
lib_raevol_X11_debug_a_SOURCES = $(common_sources) $(X11_regulation_sources)
lib_raevol_X11_debug_a_CPPFLAGS = $(AM_CPPFLAGS)

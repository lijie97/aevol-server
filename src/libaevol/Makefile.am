############################################
#                 Variables                #
############################################
SUBDIRS = SFMT-src-1.4 . tests

########################
# Set library to build #
########################
if WITH_RAEVOL
  noinst_LIBRARIES = libraevol.a
else # WITHOUT_RAEVOL
  noinst_LIBRARIES = libaevol.a
endif # RAEVOL


AM_CPPFLAGS = -I./

############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS



common_headers = ae_align.h ae_codon.h ae_dna.h ae_dump.h ae_enums.h ae_environment.h ae_env_segment.h ae_exp_manager.h ae_exp_setup.h ae_fuzzy_set.h ae_gaussian.h ae_genetic_unit.h ae_grid_cell.h ae_individual.h ae_list.h ae_logs.h ae_macros.h ae_mutation.h ae_object.h ae_output_manager.h ae_params_init.h ae_params_mut.h ae_params_record.h ae_phenotype.h ae_population.h ae_point_1d.h ae_point_2d.h ae_protein.h ae_jumping_mt.h ae_jump_poly.h ae_replication_report.h ae_dna_replic_report.h ae_rna.h ae_stats.h ae_stat_record.h ae_string.h ae_tree.h ae_utils.h ae_vis_a_vis.h ae_selection.h ae_spatial_structure.h
common_headers += params.h param_loader.h f_line.h

common_sources = ae_align.cpp ae_codon.cpp ae_dna.cpp ae_dump.cpp ae_environment.cpp ae_exp_manager.cpp ae_exp_setup.cpp ae_fuzzy_set.cpp ae_genetic_unit.cpp ae_grid_cell.cpp ae_individual.cpp ae_logs.cpp ae_mutation.cpp ae_output_manager.cpp ae_params_init.cpp ae_params_mut.cpp ae_params_record.cpp ae_phenotype.cpp ae_population.cpp ae_protein.cpp ae_jumping_mt.cpp ae_replication_report.cpp ae_dna_replic_report.cpp ae_rna.cpp ae_stats.cpp ae_stat_record.cpp ae_string.cpp ae_tree.cpp ae_vis_a_vis.cpp ae_selection.cpp ae_spatial_structure.cpp
common_sources += params.cpp param_loader.cpp f_line.cpp 


X11_headers = ae_fuzzy_set_X11.h ae_individual_X11.h ae_population_X11.h ae_exp_manager_X11.h ae_X11_window.h
X11_sources = ae_fuzzy_set_X11.cpp ae_individual_X11.cpp ae_population_X11.cpp ae_exp_manager_X11.cpp ae_X11_window.cpp

raevol_headers = ae_array_short.h ae_individual_R.h ae_influence_R.h ae_matrix_double.h ae_protein_R.h ae_rna_R.h ae_selection_R.h
raevol_sources = ae_individual_R.cpp ae_influence_R.cpp ae_protein_R.cpp ae_rna_R.cpp ae_selection_R.cpp


raevol_X11_headers = $(X11_headers) $(raevol_headers) ae_individual_R_X11.h
raevol_X11_sources = $(X11_sources) $(raevol_sources) ae_individual_R_X11.cpp




# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib


if WITH_RAEVOL
  libraevol_adir = ./
  libraevol_a_CPPFLAGS = $(AM_CPPFLAGS)
else # WITHOUT_RAEVOL
  libaevol_adir = ./
  libaevol_a_CPPFLAGS = $(AM_CPPFLAGS)
endif # RAEVOL

if WITH_RAEVOL
if WITH_X
  libraevol_a_HEADERS = $(common_headers) $(raevol_X11_headers)
  libraevol_a_SOURCES = $(common_sources) $(raevol_X11_sources)
else # WITHOUT_X
  libraevol_a_HEADERS = $(common_headers) $(raevol_headers)
  libraevol_a_SOURCES = $(common_sources) $(raevol_sources)
endif # X
else # WITHOUT_RAEVOL
if WITH_X
  libaevol_a_HEADERS = $(common_headers) $(X11_headers)
  libaevol_a_SOURCES = $(common_sources) $(X11_sources)
else # WITHOUT_X
  libaevol_a_HEADERS = $(common_headers)
  libaevol_a_SOURCES = $(common_sources)
endif # X
endif # RAEVOL

# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(SFMT-src-1.4)


# ============================================================================
# Find packages
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Add PWD to include directories
# ============================================================================
include_directories(".")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_library(aevol
  ae_align.cpp
  ae_align.h
  ae_array_short.h
  ae_codon.cpp
  ae_codon.h
  ae_dna.cpp
  ae_dna.h
  ae_dump.cpp
  ae_dump.h
  ae_enums.h
  ae_env_segment.h
  ExpManager.cpp
  ExpManager.h
  ExpManager_X11.cpp
  ExpManager_X11.h
  ExpSetup.cpp
  ExpSetup.h
  ae_gaussian.h
  ae_gene_mutation.cpp
  ae_gene_mutation.h
  ae_gene_tree.cpp
  ae_gene_tree.h
  ae_gene_tree_node.cpp
  ae_gene_tree_node.h
  GeneticUnit.cpp
  GeneticUnit.h
  habitat.cpp
  habitat.h
  ae_grid_cell.cpp
  ae_grid_cell.h
  Individual.cpp
  Individual.h
  ae_individual_X11.cpp
  ae_individual_X11.h
  individual_factory.cpp
  individual_factory.h
  ae_jumping_mt.cpp
  ae_jumping_mt.h
  ae_jump_poly.h
  ae_logs.cpp
  ae_logs.h
  ae_mutation.cpp
  ae_mutation.h
  ae_output_manager.cpp
  ae_output_manager.h
  ae_params_mut.cpp
  ae_params_mut.h
  phenotype.h
  ae_protein.cpp
  ae_protein.h
  ae_rand_mt.cpp
  ae_rand_mt.h
  ae_replication_report.cpp
  ae_replication_report.h
  ae_rna.cpp
  ae_rna.h
  ae_selection.cpp
  ae_selection.h
  ae_stat_record.cpp
  ae_stat_record.h
  ae_stats.cpp
  ae_stats.h
  ae_string.cpp
  ae_string.h
  ae_tree.cpp
  ae_tree.h
  ae_utils.h
  ae_vis_a_vis.cpp
  ae_vis_a_vis.h
  ae_X11_window.cpp
  ae_X11_window.h
  dna_replic_report.cpp
  dna_replic_report.h
  f_line.cpp
  f_line.h
  fuzzy.cpp
  fuzzy.h
  indiv_stats.h
  macros.h
  non_coding_stats.h
  param_loader.cpp
  param_loader.h
  phenotypic_target_handler.cpp
  phenotypic_target_handler.h
  phenotypic_target.cpp
  phenotypic_target.h
  point.cpp
  point.h
  Time.cpp
  Time.h
  world.cpp
  world.h)


# ============================================================================
# Use GSL and zlib
# ============================================================================
target_link_libraries(aevol PUBLIC ${ZLIB_LIBRARY})
target_link_libraries(aevol PUBLIC sfmt)


# ============================================================================
# Make STDC MACROS available (for fixed width integers)
# ============================================================================
target_compile_definitions(aevol
  PUBLIC __STDC_FORMAT_MACROS
  PUBLIC __STDC_CONSTANT_MACROS)


# ============================================================================
# Use C++11
# ============================================================================
target_compile_options(aevol PRIVATE "-std=c++14")

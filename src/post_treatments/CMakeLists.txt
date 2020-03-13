# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================


# ============================================================================
# Find packages
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Add PWD to include directories
# ============================================================================
include_directories(".")
include_directories("../libaevol")
include_directories("../libaevol/mutation")
include_directories("../libaevol/phenotype")
include_directories("../libaevol/biochemistry")
include_directories("../libaevol/rng")
include_directories("../libaevol/world")
include_directories("../libaevol/population")
include_directories("../libaevol/stats")
include_directories("../libaevol/cuda")
include_directories("../libaevol/io")
include_directories("../libaevol/gui")
include_directories("../libaevol/utils")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(aevol_post_create_eps create_eps.cpp)
target_link_libraries(aevol_post_create_eps PUBLIC aevol)

add_executable(aevol_post_extract extract.cpp)
target_link_libraries(aevol_post_extract PUBLIC aevol)

add_executable(aevol_post_ancestor_stats ancestor_stats.cpp)
target_link_libraries(aevol_post_ancestor_stats PUBLIC aevol)

add_executable(aevol_post_ancestor_mutagenesis ancestor_mutagenesis.cpp IndivAnalysis.cpp)
target_link_libraries(aevol_post_ancestor_mutagenesis PUBLIC aevol)

add_executable(aevol_post_coalescence coalescence.cpp)
target_link_libraries(aevol_post_coalescence PUBLIC aevol)

add_executable(aevol_post_protein_map protein_map.cpp)
target_link_libraries(aevol_post_protein_map PUBLIC aevol)

add_executable(aevol_post_ancestor_robustness ancestor_robustness.cpp IndivAnalysis.cpp)
target_link_libraries(aevol_post_ancestor_robustness PUBLIC aevol)

add_executable(aevol_post_lineage lineage.cpp)
target_link_libraries(aevol_post_lineage PUBLIC aevol)

add_executable(aevol_post_mutagenesis mutagenesis.cpp)
target_link_libraries(aevol_post_mutagenesis PUBLIC aevol)

add_executable(aevol_post_soft_lineage soft_lineage.cpp)
target_link_libraries(aevol_post_soft_lineage PUBLIC aevol)

add_executable(aevol_post_robustness robustness.cpp IndivAnalysis.cpp)
target_link_libraries(aevol_post_robustness PUBLIC aevol)


# ============================================================================
# Set behaviour on make install
# ============================================================================

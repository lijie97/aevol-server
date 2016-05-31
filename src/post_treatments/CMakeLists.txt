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

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(create_eps create_eps.cpp)
target_link_libraries(create_eps PUBLIC aevol)

add_executable(extract extract.cpp)
target_link_libraries(extract PUBLIC aevol)

add_executable(ancestor_stats ancestor_stats.cpp)
target_link_libraries(ancestor_stats PUBLIC aevol)

add_executable(ancestor_robustness ancestor_robustness.cpp IndivAnalysis.cpp)
target_link_libraries(mutationalrobustness PUBLIC aevol)

add_executable(lineage lineage.cpp)
target_link_libraries(lineage PUBLIC aevol)

add_executable(mutagenesis mutagenesis.cpp)
target_link_libraries(mutagenesis PUBLIC aevol)


# ============================================================================
# Set behaviour on make install
# ============================================================================

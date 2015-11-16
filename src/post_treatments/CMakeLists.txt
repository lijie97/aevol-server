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

add_executable(ancstats ancstats.cpp)
target_link_libraries(ancstats PUBLIC aevol)

add_executable(lineage lineage.cpp)
target_link_libraries(lineage PUBLIC aevol)

add_executable(fixed_mutations fixed_mutations.cpp)
target_link_libraries(fixed_mutations PUBLIC aevol)

add_executable(mutagenesis mutagenesis.cpp)
target_link_libraries(mutagenesis PUBLIC aevol)

add_executable(mutationalrobustness mutationalrobustness.cpp)
target_link_libraries(mutationalrobustness PUBLIC aevol)

add_executable(robustness robustness.cpp)
target_link_libraries(robustness PUBLIC aevol)

add_executable(transform_plasmid transform_plasmid.cpp)
target_link_libraries(transform_plasmid PUBLIC aevol)

add_executable(view_generation view_generation.cpp)
target_link_libraries(view_generation PUBLIC aevol)

# ============================================================================
# Set behaviour on make install
# ============================================================================

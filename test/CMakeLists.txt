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

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(TranscriptionTranslationTest TranscriptionTranslationTest.cpp)
target_link_libraries(TranscriptionTranslationTest PUBLIC aevol)

add_executable(IndividualTest IndividualTest.cpp)
target_link_libraries(IndividualTest PUBLIC aevol)


# ============================================================================
# Set behaviour on make install
# ============================================================================

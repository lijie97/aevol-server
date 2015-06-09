# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(libaevol)


# ============================================================================
# Find packages
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Add PWD to include directories
# ============================================================================
include_directories(".")
include_directories("libaevol")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(aevol_create aevol_create.cpp)
add_executable(aevol_run aevol_run.cpp)

target_link_libraries(aevol_create PUBLIC aevol)
target_link_libraries(aevol_create PUBLIC sfmt)
target_link_libraries(aevol_run PUBLIC aevol)
target_link_libraries(aevol_run PUBLIC sfmt)

# ============================================================================
# Use zlib
# ============================================================================
target_link_libraries(aevol_create PUBLIC ${ZLIB_LIBRARY})
target_link_libraries(aevol_run PUBLIC ${ZLIB_LIBRARY})


# ============================================================================
# Make STDC MACROS available (for fixed width integers)
# ============================================================================
target_compile_definitions(aevol_create
  PUBLIC __STDC_FORMAT_MACROS
  PUBLIC __STDC_CONSTANT_MACROS)
target_compile_definitions(aevol_run
  PUBLIC __STDC_FORMAT_MACROS
  PUBLIC __STDC_CONSTANT_MACROS)


# ============================================================================
# Use C++11
# ============================================================================
target_compile_options(aevol_create PRIVATE "-std=c++14")
target_compile_options(aevol_run PRIVATE "-std=c++14")


# ============================================================================
# Set behaviour on make install
# ============================================================================
install(TARGETS aevol_create DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_create)
install(TARGETS aevol_run DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_run)


# ============================================================================
# Include tests
# ============================================================================
#add_subdirectory(tests EXCLUDE_FROM_ALL)

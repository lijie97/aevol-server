# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(libaevol)
add_subdirectory(post_treatments)


# ============================================================================
# Find packages
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Add PWD to include directories
# ============================================================================
include_directories(".")
include_directories("libaevol")
include_directories("libaevol/mutation")
include_directories("libaevol/phenotype")
include_directories("libaevol/biochemistry")
include_directories("libaevol/rng")
include_directories("libaevol/world")
include_directories("libaevol/population")
include_directories("libaevol/stats")
include_directories("libaevol/cuda")
include_directories("libaevol/io")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(aevol_create aevol_create.cpp)
target_link_libraries(aevol_create PUBLIC aevol)

add_executable(aevol_modify aevol_modify.cpp)
target_link_libraries(aevol_modify PUBLIC aevol)

#######################
## Searching for libs
#######################
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.65.0 COMPONENTS filesystem)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(aevol_propagate aevol_propagate.cpp)
    target_link_libraries(aevol_propagate PUBLIC aevol ${Boost_LIBRARIES})
endif()


add_executable(aevol_run aevol_run.cpp)
target_link_libraries(aevol_run PUBLIC aevol)


# ============================================================================
# Set behaviour on make install
# ============================================================================
install(TARGETS aevol_create DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_create)
install(TARGETS aevol_modify DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_modify)
install(TARGETS aevol_propagate DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_propagate)
install(TARGETS aevol_run DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_run)


# ============================================================================
# Include tests
# ============================================================================
add_subdirectory(tests EXCLUDE_FROM_ALL)

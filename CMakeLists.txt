# ============================================================================
# Require minimal version of cmake
# ============================================================================
cmake_minimum_required(VERSION 3.0.2)


# ============================================================================
# Set project name and languages
# ============================================================================
project(aevol VERSION 5.0 LANGUAGES C CXX)
add_definitions(-DVERSION="${PROJECT_VERSION}")

# Temporary stuff to be fully integrated into CMake
add_definitions(-DSFMT_MEXP=607)
add_definitions(-D__X11)
add_definitions(-DDEBUG)

FIND_PACKAGE ( X11 REQUIRED )
IF ( X11_FOUND )
      INCLUDE_DIRECTORIES ( ${X11_INCLUDE_DIR} )
      LINK_LIBRARIES ( ${X11_LIBRARIES} )
ENDIF ( X11_FOUND )


# ============================================================================
# Get GNU standard installation directories (GNUInstallDirs module)
# ============================================================================
include(GNUInstallDirs)


# ============================================================================
# Tell CMake where to look for custom modules
# ============================================================================
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# ============================================================================
# Tell cmake where to put binary files.
# By GNU standards "executable programs that users can run" should go in
# bindir a.k.a ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
# and "executable programs to be run by other programs rather than by users"
# in libexecdir a.k.a ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}
# ============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# ============================================================================
# Set build type specific compilation flags
# ============================================================================
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")


# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(src)


# ============================================================================
# Adds the 'dist' target (that will use CPack)
# ============================================================================
#add_custom_target(dist COMMAND ${CMAKE_BUILD_TOOL} package_source)


# ============================================================================
# Add the 'uninstall' target (uses a custom script)
# ============================================================================
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake)

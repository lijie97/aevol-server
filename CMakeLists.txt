# ============================================================================
# Require minimal version of cmake
# ============================================================================
cmake_minimum_required(VERSION 3.0.2)


# ============================================================================
# Set project name and languages
# ============================================================================
project(aevol VERSION 5.0 LANGUAGES C CXX)
add_definitions(-DVERSION="${PROJECT_VERSION}")

set(AUTHORIZED_MTPERIOD 607 1279 2281 4253 11213 19937 44497 86243 132049 216091)


set(with-x OFF CACHE BOOL "Whether to enable graphical outputs")
set(with-omp OFF CACHE BOOL "Whether to enable OpenMP parallelization")
set(with-tbb OFF CACHE BOOL "Whether to enable TBB parallelization")
set(enable-profiling OFF CACHE BOOL "Whether to enable profiling")
set(with-tracing OFF CACHE BOOL "Whether to use tracing")
set(enable-perflog OFF CACHE BOOL "Whether to active performance log")
set(enable-normalized-fitness OFF CACHE BOOL "With this option, the NORMALIZED_FITNESS flag is defined, allowing a different fitness calculation")
set(enable-mtperiod 607 CACHE STRING "Period of the Mersen Twister. Autorized values are : 607 1279 2281 4253 11213 19937 44497 86243 132049 216091")
set(enable-trivialjumps OFF CACHE STRING "When this option is set, a trivial jump algorithm will be used instead of the polynomial-based method")

set(with-triangle OFF CACHE BOOL "Whether to enable triangle phenotypic target (else Gaussian)")
set(with-raevol OFF CACHE BOOL "Whether to enable R-Aevol (Gene regulation) version")
set(with-detectclone ON CACHE BOOL "Whether to enable clones and not recompute them")
set(with-floatconcentration OFF CACHE BOOL "Whether to enable the encoding of concentration has float (and not double)")
set(with-simd OFF CACHE BOOL "Whether to switch to the SIMD implementation of (R-)Aevol")
set(with-perf-traces OFF CACHE BOOL "Whether to activate performance traces of (R-)Aevol")
set(search-type "vanilla" CACHE STRING "Which SIMD pattern search version of (R-)Aevol")

if ( ${with-x} )
    FIND_PACKAGE(X11 REQUIRED)
    if ( X11_FOUND )
        add_definitions(-D__X11)
        INCLUDE_DIRECTORIES(${X11_INCLUDE_DIR})
        LINK_LIBRARIES(${X11_LIBRARIES})
        message("graphical output enable")
    endif ( X11_FOUND )
else ()
    add_definitions(-D__NO_X)
    message("graphical output disable")
endif ()

if ( ${with-omp} )
    find_package(OpenMP REQUIRED)
    if ( OPENMP_FOUND )
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()

if ( ${enable-profiling} )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif ()

if ( ${with-tracing} )
    add_definitions(-D__TRACING__)
endif ()

if ( ${enable-perflog} )
    add_definitions(-D__PERF_LOG__)
endif ()

if ( ${search-type} STREQUAL "vanilla")
    add_definitions(-DVANILLA_SEARCH)
elseif ( ${search-type} STREQUAL "cmod")
    add_definitions(-DCMOD_SEARCH)
elseif ( ${search-type} STREQUAL "amod")
    add_definitions(-DAMOD_SEARCH)
elseif ( ${search-type} STREQUAL "b1mod")
    add_definitions(-DB1MOD_SEARCH)
elseif ( ${search-type} STREQUAL "b2mod")
    add_definitions(-DB2MOD_SEARCH) 
endif()
message("Search type is set to ${search-type}")


if ( ${enable-normalized-fitness} )
    add_definitions(-DNORMALIZED_FITNESS)
endif ()

if ( DEFINED enable-mtperiod )
    list(FIND AUTHORIZED_MTPERIOD ${enable-mtperiod} _index)
    if ( ${_index} GREATER -1 )
        set(MTPERIOD ${enable-mtperiod})
    else ()
        message(FATAL_ERROR "period is not a valid Mersenne Twister period")
    endif ()
endif ()
message("Mersen Twister period is set to ${MTPERIOD}")
add_definitions(-DSFMT_MEXP=${MTPERIOD})

if ( ${enable-trivialjumps} )
    if ( ${enable-trivialjumps} MATCHES "[0-9]+" )
        add_definitions(-DTRIVIAL_METHOD_JUMP_SIZE=${enable-trivialjumps})
    else ()
        add_definitions(-DTRIVIAL_METHOD_JUMP_SIZE=1000)
    endif ()
endif ()

if ( ${with-triangle} )
    add_definitions(-DPHENOTYPIC_TARGET_TRIANGLE)
endif ()

if ( ${with-raevol} )
    add_definitions(-D__REGUL)
endif ()

if ( ${with-detectclone} )
    add_definitions(-D__DETECT_CLONE)
endif ()

if ( ${with-floatconcentration} )
    add_definitions(-D__FLOAT_CONCENTRATION)
endif ()

if ( ${with-simd} )
    add_definitions(-DWITH_STANDALONE_SIMD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=skylake -O3 -AVX2 -g")
endif ()

if ( ${with-perf-traces} )
    add_definitions(-DWITH_PERF_TRACES)
endif ()
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
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -Wall -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -Wall -DDEBUG")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -ffast-math -pipe -fexpensive-optimizations")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -AVX512 -march=skylake -O3 -ffast-math -pipe -fexpensive-optimizations")


# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(src)

if(APPLE)
   include_directories(AFTER "/opt/X11/include")
endif()

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

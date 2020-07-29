# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(libaevol)

if (with-post-treatments)
  add_subdirectory(post_treatments)
endif()


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
include_directories("libaevol/gui")
include_directories("libaevol/utils")

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

if ( ${with-simd} )
    add_definitions(-DWITH_STANDALONE_SIMD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()

if ( USE_CUDA )
    find_package(CUDA REQUIRED)
    enable_language(CUDA)
    set(CUDA_NVCC_FLAGS "-std=c++11 -arch=compute_60 -lineinfo")
    cuda_compile(cuda_aevol libaevol/cuda/Algorithms.cu libaevol/cuda/ThreefryGPU.cu)
endif ( USE_CUDA )

if ( CUDA_FOUND )
    add_definitions(-DUSE_CUDA)
    cuda_add_executable(aevol_run_gpu aevol_run.cpp ${cuda_aevol})
    # nvToolsExt for enhanced profiling (ad-hoc chunks)
    target_link_libraries(aevol_run_gpu aevol nvToolsExt)
    install(TARGETS aevol_run_gpu DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_run_gpu)
else ( CUDA_FOUND )
    add_executable(aevol_run aevol_run.cpp)
    target_link_libraries(aevol_run PUBLIC aevol)
    install(TARGETS aevol_run DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_run)
endif ( CUDA_FOUND )

# ============================================================================
# Set behaviour on make install
# ============================================================================
install(TARGETS aevol_create DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_create)
install(TARGETS aevol_modify DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_modify)
install(TARGETS aevol_propagate DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_propagate)

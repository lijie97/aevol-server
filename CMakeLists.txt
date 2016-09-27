cmake_minimum_required(VERSION 3.3)
project(raevol_cuda)

set(__OPENMP_GPU)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp src/libaevol/ProxyIndividual.cpp src/libaevol/ProxyIndividual.h src/libaevol/algorithm_cuda.cpp src/libaevol/algorithm_cuda.h src/generate_lookup_map_pow.cpp)
add_executable(raevol_cuda ${SOURCE_FILES})

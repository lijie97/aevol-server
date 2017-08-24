cmake_minimum_required(VERSION 3.3)
project(raevol_cuda)

set(__OPENMP_GPU)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp src/libaevol/ProxyIndividual.cpp src/libaevol/ProxyIndividual.h src/libaevol/algorithm_cuda.cpp src/libaevol/algorithm_cuda.h src/generate_lookup_map_pow.cpp src/regulation_network.cpp src/ode_solver_only.cpp src/ode_solver_only.h src/ode_solver_only_gpu.h src/ode_solver_only_omp_gpu.cpp src/ode_solver_only_omp_gpu.h src/rk_solver.cpp src/rk_solver.h src/libaevol/CUDA_Individual_cu.h src/libaevol/SIMD_Individual.cpp src/libaevol/SIMD_Individual.h src/libaevol/Dna_SIMD.cpp src/libaevol/Dna_SIMD.h src/libaevol/MutationEvent.cpp src/libaevol/MutationEvent.h src/libaevol/MutatePromoterEvent.cpp src/libaevol/MutatePromoterEvent.h)
add_executable(raevol_cuda ${SOURCE_FILES})

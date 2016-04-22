cmake_minimum_required(VERSION 3.3)
project(raevol_cuda)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp src/libaevol/ProxyIndividual.cpp src/libaevol/ProxyIndividual.h)
add_executable(raevol_cuda ${SOURCE_FILES})

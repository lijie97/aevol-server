cmake_minimum_required(VERSION 3.2)
project(aevol_regul)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(aevol_regul ${SOURCE_FILES})
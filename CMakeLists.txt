cmake_minimum_required(VERSION 2.8.12)
project(EM_C++)

include_directories(include/)

set(CMAKE_CXX_FLAGS "-O3")

add_executable(fastppca src/fastppca.cpp src/genotype.cpp src/storage.cpp)

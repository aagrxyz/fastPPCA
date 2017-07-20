cmake_minimum_required(VERSION 2.8.12)
project(EM_C++)

include_directories(include)

set(CMAKE_CXX_FLAGS "-O3 -g -march=native")

# add_executable(em_naive src/em_naive.cpp src/genotype.cpp src/storage.cpp)
# add_executable(fastppca src/em_mailman.cpp src/genotype.cpp src/storage.cpp)
add_executable(fastppca src/fastppca.cpp src/genotype.cpp src/storage.cpp)

cmake_minimum_required(VERSION 2.8.12)
project(EM_C++)

# set(EIGEN_INCLUDE_DIR "" CACHE FILEPATH "Path to Eigen includes")

# include_directories(${EIGEN_INCLUDE_DIR})
include_directories(include)

set(CMAKE_CXX_FLAGS "-O3 -march=native")

add_executable(em_naive src/em_naive.cpp src/genotype.cpp)
add_executable(em_mailman src/em_mailman.cpp src/genotype.cpp)
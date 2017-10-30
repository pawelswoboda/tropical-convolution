project(tropical-convolution)

cmake_minimum_required(VERSION 2.8.12)

# C++14
add_compile_options(-std=c++14)

add_executable(tropical_convolution_test tropical_convolution_test.cpp)
add_test(NAME TROPICAL_CONVOLUTION_TEST COMMAND tropical_convolution_test)

add_executable(tropical_convolution_benchmark tropical_convolution_benchmark.cpp)


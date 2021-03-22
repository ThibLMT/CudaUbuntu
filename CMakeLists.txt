cmake_minimum_required(VERSION 3.17)
project(CudaUbuntu CUDA)

set(CMAKE_CUDA_STANDARD 14)

add_executable(CudaUbuntu main.cu def_types.h def_const.h def_global_variables.h)

set_target_properties(
        CudaUbuntu
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
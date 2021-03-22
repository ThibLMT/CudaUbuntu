cmake_minimum_required(VERSION 3.17)
project(CudaUbuntu CUDA)

set(CMAKE_CUDA_STANDARD 14)

add_executable(CudaUbuntu main.cu def_type.h def_const.h)

set_target_properties(
        CudaUbuntu
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
cmake_minimum_required(VERSION 3.17)
project(CudaUbuntu CUDA)

set(CMAKE_CUDA_STANDARD 14)

add_executable(CudaUbuntu main.cu def_types.h def_const.h def_global_variables.h allocate_variables.cu allocate_variables.cuh
        read_micro.cuh read_micro.cu init_params.cu init_params.cuh sub_domain.cu sub_domain.cuh write_micro.cu write_micro.cuh
        cal_interactions_forces.cu cal_interactions_forces.cuh particle_interactions.cu particle_interactions.cuh
        Lock.cuh helper_cuda.h helper_functions.h tangential_force_mindlin_vect.cu tangential_force_mindlin_vect.cuh small_vect_rot.cu small_vect_rot.cuh rolling_resistance_torque.cu rolling_resistance_torque.cuh normal_force_hertz.cu normal_force_hertz.cuh)

set_target_properties(
        CudaUbuntu
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
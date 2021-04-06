//
// Created by ThibLMT on 06/04/2021.
//

#ifndef CUDAUBUNTU_LOCK_CUH
#define CUDAUBUNTU_LOCK_CUH


struct Lock {

    int *d_state;

    // --- Constructor
    Lock(void) {
        int h_state = 0;                                        // --- Host side lock state initializer
        cudaMalloc((void **)&d_state, sizeof(int));  // --- Allocate device side lock state
        cudaMemcpy(d_state, &h_state, sizeof(int), cudaMemcpyHostToDevice); // --- Initialize device side lock state
    }

    // --- Destructor
    __host__ __device__ ~Lock(void) {
        cudaFree(d_state);
    }

    // --- Lock function
    __device__ void lock(void) { while (atomicCAS(d_state, 0, 1) != 0); }

    // --- Unlock function
    __device__ void unlock(void) { atomicExch(d_state, 0); }
};

#endif //CUDAUBUNTU_LOCK_CUH

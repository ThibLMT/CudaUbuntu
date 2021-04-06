//
// Created by ThibLMT on 06/04/2021.
//

#ifndef CUDAUBUNTU_LOCK_CUH
#define CUDAUBUNTU_LOCK_CUH

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

struct Lock {

    int *d_state;

    // --- Constructor
    Lock(void) {
        int h_state = 0;                                        // --- Host side lock state initializer
        gpuErrchk(cudaMalloc((void **)&d_state, sizeof(int)));  // --- Allocate device side lock state
        gpuErrchk(cudaMemcpy(d_state, &h_state, sizeof(int), cudaMemcpyHostToDevice)); // --- Initialize device side lock state
    }

    // --- Destructor
    __host__ __device__ ~Lock(void) {
        gpuErrchk(cudaFree(d_state));
    }

    // --- Lock function
    __device__ void lock(void) { while (atomicCAS(d_state, 0, 1) != 0); }

    // --- Unlock function
    __device__ void unlock(void) { atomicExch(d_state, 0); }
};

#endif //CUDAUBUNTU_LOCK_CUH

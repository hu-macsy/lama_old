#pragma once

#ifdef __CUDACC__
    #define CUDA_CALLABLE_MEMBER __device__ __host__
#else
    #define CUDA_CALLABLE_MEMBER
#endif

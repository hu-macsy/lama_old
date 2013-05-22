/**
 * @file CUDAError.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class CUDAError.
 * @author Thomas Brandes
 * @date 15.07.2011
 * $Id$
 */

// hpp
#include <lama/cuda/CUDAError.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>

namespace lama
{

const char* cudaDriverErrorString( CUresult res )
{
    const char* str = "";

    switch ( res )

    {
    case CUDA_SUCCESS:
        str = "No errors";
        break;
    case CUDA_ERROR_INVALID_VALUE:
        str = "Invalid value.";
        break;
    case CUDA_ERROR_OUT_OF_MEMORY:
        str = "Out of memory.";
        break;
    case CUDA_ERROR_NOT_INITIALIZED:
        str = "Driver not initialized.";
        break;
    case CUDA_ERROR_DEINITIALIZED:
        str = "Driver deinitialized.";
        break;
    case CUDA_ERROR_NO_DEVICE:
        str = "No CUDA-capable device available.";
        break;
    case CUDA_ERROR_INVALID_DEVICE:
        str = "Invalid device.";
        break;
    case CUDA_ERROR_INVALID_IMAGE:
        str = "Invalid kernel image.";
        break;
    case CUDA_ERROR_INVALID_CONTEXT:
        str = "Invalid context.";
        break;
    case CUDA_ERROR_CONTEXT_ALREADY_CURRENT:
        str = "Context already current.";
        break;
    case CUDA_ERROR_MAP_FAILED:
        str = "Map failed.";
        break;
    case CUDA_ERROR_UNMAP_FAILED:
        str = "Unmap failed.";
        break;
    case CUDA_ERROR_ARRAY_IS_MAPPED:
        str = "Array is mapped.";
        break;
    case CUDA_ERROR_ALREADY_MAPPED:
        str = "Already mapped.";
        break;
    case CUDA_ERROR_NO_BINARY_FOR_GPU:
        str = "No binary for GPU.";
        break;
    case CUDA_ERROR_ALREADY_ACQUIRED:
        str = "Already acquired.";
        break;
    case CUDA_ERROR_NOT_MAPPED:
        str = "Not mapped.";
        break;
    case CUDA_ERROR_INVALID_SOURCE:
        str = "Invalid source.";
        break;
    case CUDA_ERROR_FILE_NOT_FOUND:
        str = "File not found.";
        break;
    case CUDA_ERROR_INVALID_HANDLE:
        str = "Invalid handle.";
        break;
    case CUDA_ERROR_NOT_FOUND:
        str = "Not found.";
        break;
    case CUDA_ERROR_NOT_READY:
        str = "CUDA not ready.";
        break;
    case CUDA_ERROR_LAUNCH_FAILED:
        str = "Launch failed.";
        break;
    case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES:
        str = "Launch exceeded resources.";
        break;
    case CUDA_ERROR_LAUNCH_TIMEOUT:
        str = "Launch exceeded timeout.";
        break;
    case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING:
        str = "Launch with incompatible texturing.";
        break;
    case CUDA_ERROR_UNKNOWN:
        str = "Unknown error. ";
        break;
    default:
        str = "Illegal result value error.";
    }

    return str;
}

const char* cublasErrorString( cublasStatus res )
{
    const char* str = "";

    switch ( res )

    {
    case CUBLAS_STATUS_SUCCESS:
        str = "CUBLAS successful";
        break;
    case CUBLAS_STATUS_NOT_INITIALIZED:
        str = "CUBLAS library not initialized";
        break;
    case CUBLAS_STATUS_ALLOC_FAILED:
        str = "resource allocation failed";
        break;
    case CUBLAS_STATUS_INVALID_VALUE:
        str = "unsupported numerical value was passed to function";
        break;
    case CUBLAS_STATUS_ARCH_MISMATCH:
        str = "function requires an architectural feature absent from the architecture of the device";
        break;
    case CUBLAS_STATUS_MAPPING_ERROR:
        str = "access to GPU memory space failed";
        break;
    case CUBLAS_STATUS_EXECUTION_FAILED:
        str = "GPU program failed to execute";
        break;
    case CUBLAS_STATUS_INTERNAL_ERROR:
        str = "CUBLAS library not initialized";
        break;
    default:
        str = "Unknown CUBLAS error";
    }

    return str;
}

/* ----------------------------------------------------------------------------- */

} //namespace lama


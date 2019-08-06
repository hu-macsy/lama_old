/**
 * @file CUDAError.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains the implementation of the class CUDAError.
 * @author Thomas Brandes
 * @date 15.07.2011
 */

// hpp
#include <scai/common/cuda/CUDAError.hpp>

// CUDA
#include <cuda.h>
#if (defined(__clang__) && __clang_major__ == 3 && __clang_minor__ <= 3)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wunneeded-internal-declaration"
#include <cuda_runtime.h>
# pragma GCC diagnostic pop
#else
#include <cuda_runtime.h>
#endif
#include <cusparse.h>

namespace scai
{

namespace common
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

        case CUDA_ERROR_ILLEGAL_ADDRESS:
            str = "Illegal Address (inconsistent state, process should terminate)";
            break;

        default:
            str = "Illegal result value error.";
    }

    return str;
}

const char* cublasErrorString( cublasStatus_t res )
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

const char* cusparseErrorString( cusparseStatus_t res )
{
    const char* str = "";

    switch ( res )
    {
        case CUSPARSE_STATUS_SUCCESS:
            str = "CUSPARSE successful";
            break;

        case CUSPARSE_STATUS_NOT_INITIALIZED:
            str = "CUSPARSE library not initialized";
            break;

        case CUSPARSE_STATUS_ALLOC_FAILED:
            str = "resource allocation failed";
            break;

        case CUSPARSE_STATUS_INVALID_VALUE:
            str = "unsupported numerical value was passed to function";
            break;

        case CUSPARSE_STATUS_ARCH_MISMATCH:
            str = "function requires an architectural feature absent from the architecture of the device";
            break;

        case CUSPARSE_STATUS_MAPPING_ERROR:
            str = "access to GPU memory space failed";
            break;

        case CUSPARSE_STATUS_EXECUTION_FAILED:
            str = "GPU program failed to execute";
            break;

        case CUSPARSE_STATUS_INTERNAL_ERROR:
            str = "CUSPARSE internal error";
            break;

        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            str = "CUSPARSE unsupported matrix type";
            break;

        case CUSPARSE_STATUS_ZERO_PIVOT:
            str = "CUSPARSE zero pivot";
            break;

        default:
            str = "Unknown CUSPARSE error";
    }

    return str;
}

#if ( CUDART_VERSION >= 7050 )
const char* cusolverErrorString( cusolverStatus_t res )
{
    const char* str = "";

    switch ( res )
    {
        case CUSOLVER_STATUS_SUCCESS:
            str = "CUSOLVER successful";
            break;

        case CUSOLVER_STATUS_NOT_INITIALIZED:
            str = "CUSOLVER library not initialized";
            break;
        case CUSOLVER_STATUS_ALLOC_FAILED:
            str = "resource allocation failed";
            break;
        case CUSOLVER_STATUS_INVALID_VALUE:
            str = "unsupported numerical value was passed to function";
            break;
        case CUSOLVER_STATUS_ARCH_MISMATCH:
            str = "function requires an architectural feature absent from the architecture of the device";
            break;
        case CUSOLVER_STATUS_EXECUTION_FAILED:
            str = "CUSOLVER program failed to execute";
            break;
        case CUSOLVER_STATUS_INTERNAL_ERROR:
            str = "CUSOLVER internal error";
            break;
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            str = "CUSOLVER unsupported matrix type";
            break;

        default:
            str = "Unknown CUSPARSE error";
    }

    return str;
}

#endif

const char* cufftErrorString( cufftResult res )
{
    const char* str = "";

    switch ( res )
    {
        case CUFFT_SUCCESS:
            str = "CUFFT successful";
            break;
        case CUFFT_INVALID_PLAN:
            str = "CUFFT invalid plan";
            break;
        case CUFFT_ALLOC_FAILED:
            str = "CUFFT alloc failed";
            break;
        case CUFFT_INVALID_TYPE:
            str = "CUFFT invalid type";
            break;
        case CUFFT_INVALID_VALUE:
            str = "CUFFT invalid value";
            break;
        case CUFFT_INTERNAL_ERROR:
            str = "CUFFT internal error";
            break;
        case CUFFT_EXEC_FAILED:
            str = "CUFFT execution failed";
            break;
        case CUFFT_SETUP_FAILED:
            str = "CUFFT setup failed";
            break;
        case CUFFT_INVALID_SIZE:
            str = "CUFFT invalid size";
            break;
        case CUFFT_UNALIGNED_DATA:
            str = "CUFFT unaligned data";
            break;
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            str = "CUFFT incomplete parameter list";
            break;
        case CUFFT_INVALID_DEVICE:
            str = "CUFFT invalid device";
            break;
        case CUFFT_PARSE_ERROR:
            str = "CUFFT parse error";
            break;
        case CUFFT_NO_WORKSPACE:
            str = "CUFFT no workspace";
            break;
        case CUFFT_NOT_IMPLEMENTED:
            str = "CUFFT not implemented";
            break;
        case CUFFT_LICENSE_ERROR:
            str = "CUFFT license error";
            break;

        default:
            str = "Unknown CUFFT error";
        }

    return str;
}

/* ----------------------------------------------------------------------------- */

} /* end namespace common */

} /* end namespace scai */


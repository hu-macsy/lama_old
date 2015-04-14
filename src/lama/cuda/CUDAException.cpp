/**
 * @file CUDAException.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief CUDAException.cpp
 * @author Jiri Kraus
 * @date 20.05.2011
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDAException.hpp>

#include <stdio.h>
#include <sstream>

#if __GNUC__ >= 4 &&  __GNUC_MINOR__ > 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif //__GNUC__ >= 4 &&  __GNUC_MINOR__ > 6
#include <cuda_runtime.h> /* no diagnostic for this one */
#if __GNUC__ >= 4 &&  __GNUC_MINOR__ > 6
#pragma GCC diagnostic pop
#endif //__GNUC__ >= 4 &&  __GNUC_MINOR__ > 6
#include <cublas_v2.h>

namespace lama
{

CUDAException::CUDAException( const std::string& message, const cudaError_t cudaError )
{
    std::ostringstream oss;
    oss << message << " Cause: " << cudaError;
    oss << " (" << cudaGetErrorString( cudaError ) << ")";
    mMessage = oss.str();

    LAMA_LOG_WARN( logger, "EXCEPTION: " << message )
}

CUDAException::~CUDAException() throw ()
{
}

}

void lama_printError_cuda( const int error, const char* file, const int line )
{
    if( error != cudaSuccess )
    {
        fprintf( stderr, "CUDA error: %s (%d) in %s on line %d\n", cudaGetErrorString( (cudaError_t) error ), error,
                 file, line );
    }
}

void lama_printError_cublas( const int error, const char* file, const int line )
{
    char errorMsg[256];

    if( error != CUBLAS_STATUS_SUCCESS )
    {
        switch( error )
        {
            case CUBLAS_STATUS_NOT_INITIALIZED:
                sprintf( errorMsg, "CUBLAS library not initialized" );
                break;

            case CUBLAS_STATUS_ALLOC_FAILED:
                sprintf( errorMsg, "resource allocation failed" );
                break;

            case CUBLAS_STATUS_INVALID_VALUE:
                sprintf( errorMsg, "unsupported numerical value was passed to function" );
                break;

            case CUBLAS_STATUS_ARCH_MISMATCH:
                sprintf( errorMsg,
                         "function requires an architectural feature absent from the architecture of the device" );
                break;

            case CUBLAS_STATUS_MAPPING_ERROR:
                sprintf( errorMsg, "access to GPU memory space failed" );
                break;

            case CUBLAS_STATUS_EXECUTION_FAILED:
                sprintf( errorMsg, "GPU program failed to execute" );
                break;

            case CUBLAS_STATUS_INTERNAL_ERROR:
                sprintf( errorMsg, "CUBLAS library not initialized" );
                break;

            default:
                sprintf( errorMsg, "Unknown" );
        }

        fprintf( stderr, "CUBLAS error: %s (%d) in %s on line %d\n", errorMsg, error, file, line );
    }

}

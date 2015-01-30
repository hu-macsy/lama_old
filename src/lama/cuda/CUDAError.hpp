/**
 * @file CUDAError.hpp
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
 * @brief Error handling for CUDA ( runtime, api, cublas, cusparse )
 * @author Thomas Brandes
 * @date 15.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_CUDA_ERROR_HPP_
#define LAMA_CUDA_ERROR_HPP_

#include <lama/exception/Exception.hpp>
#include <lama/exception/LAMAAssert.hpp>

#include <cuda.h>
#include <cublas_v2.h>

namespace lama
{
/** Function that translates enum CUresult to strings. */

const char* cudaDriverErrorString( CUresult res );

/** Function that translates enum cublasStatus to strings. */

const char* cublasErrorString( cublasStatus_t res );

}

/** Macro for CUDA driver API calls to catch errors */

#define LAMA_CUDA_DRV_CALL(call, msg)                                               \
    {                                                                                   \
        CUresult res = call;                                                            \
        if ( CUDA_SUCCESS != res )                                                      \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "CUDA driver error in line " << __LINE__;                       \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Error: ";                                                    \
            errorStr << lama::cudaDriverErrorString( res );                             \
            errorStr << ", CUresult = " << res << "\n";                                 \
            lama::Exception::addCallStack( errorStr );                                  \
            fprintf(stderr, "%s\n", errorStr.str().c_str() );                           \
            throw lama::Exception( errorStr.str() );                                    \
        }                                                                               \
    }

#define LAMA_CUDA_RT_CALL(call, msg)                                                \
    {                                                                                   \
        cudaError_t res = call;                                                         \
        if ( cudaSuccess != res )                                                       \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "CUDA runtime error in line " << __LINE__;                      \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Error: ";                                                    \
            errorStr << cudaGetErrorString( res );                                      \
            errorStr << ", cudaError_t = " << res << "\n";                              \
            lama::Exception::addCallStack( errorStr );                                  \
            fprintf(stderr, "%s\n", errorStr.str().c_str() );                           \
            throw lama::Exception( errorStr.str() );                                    \
        }                                                                               \
    }

#define LAMA_CUBLAS_CALL( call, msg )                                               \
    {                                                                                   \
        cublasStatus_t res = call;                                                        \
        if ( CUBLAS_STATUS_SUCCESS != res )                                             \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "CUBLAS error in line " << __LINE__;                            \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Error: ";                                                    \
            errorStr << cublasErrorString( res );                                       \
            errorStr << ", cublasStatus = " << res << "\n";                             \
            lama::Exception::addCallStack( errorStr );                                  \
            fprintf(stderr, "%s\n", errorStr.str().c_str() );                           \
            throw lama::Exception( errorStr.str() );                                    \
        }                                                                               \
    }

#define LAMA_CUSPARSE_CALL( call, msg )                                             \
    {                                                                                   \
        cusparseStatus_t res = call;                                                    \
        if ( CUSPARSE_STATUS_SUCCESS != res )                                           \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "CUSparse error in line " << __LINE__;                          \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Error: ";                                                    \
            errorStr << "cusparseStatus = " << res << "\n";                             \
            lama::Exception::addCallStack( errorStr );                                  \
            fprintf(stderr, "%s\n", errorStr.str().c_str() );                           \
            throw lama::Exception( errorStr.str() );                                    \
        }                                                                               \
    }

#define LAMA_CHECK_CUDA_ACCESS                                                        \
    {                                                                                     \
        CUcontext pctx;                                                                   \
        const int cudaErrorValue = cuCtxGetCurrent( &pctx );                              \
        LAMA_ASSERT_EQUAL_ERROR( cudaErrorValue, cudaSuccess )                            \
        LAMA_ASSERT_ERROR( pctx, "No current context, forgotten LAMA_CONTEXT_ACCESS ?" ) \
    }

#define LAMA_CHECK_CUDA_ERROR                                                         \
    {                                                                                 \
        LAMA_CUDA_RT_CALL( cudaGetLastError(), "last CUDA error" )                    \
    }

//    #define LAMA_CHECK_CUDA_ERROR

#endif //  LAMA_CUDA_ERROR_HPP_

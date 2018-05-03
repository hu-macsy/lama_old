/**
 * @file CUDAError.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Error handling for CUDA ( runtime, api, cublas, cusparse )
 * @author Thomas Brandes
 * @date 15.07.2012
 */
#pragma once

// local library
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/cuda/CUDAException.hpp>

// CUDA
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <cufft.h>

#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif ( CUDART_VERSION >= 7050 )
#include <cusolverSp.h>
#endif

namespace scai
{

namespace common
{
/** Function that translates enum CUresult to strings. */

COMMON_DLL_IMPORTEXPORT const char* cudaDriverErrorString( CUresult res );

/** Function that translates enum cublasStatus to strings. */

COMMON_DLL_IMPORTEXPORT const char* cublasErrorString( cublasStatus_t res );

/** Function that translates enum cusparseStatus to strings. */

COMMON_DLL_IMPORTEXPORT const char* cusparseErrorString( cusparseStatus_t res );

/** Function that translates enum cufftResult to strings. */

COMMON_DLL_IMPORTEXPORT const char* cufftErrorString( cufftResult res );

#if ( CUDART_VERSION >= 7050 )
/** Function that translates enum cusolverStatus to strings. */

COMMON_DLL_IMPORTEXPORT const char* cusolverErrorString( cusolverStatus_t res );
#endif

} /* end namespace common */

} /* end namespace scai */

/** Macro for CUDA driver API calls to catch errors */

#define SCAI_CUDA_DRV_CALL_EXCEPTION( call, msg, exception )                        \
    {                                                                               \
        CUresult res = call;                                                        \
        if ( CUDA_SUCCESS != res )                                                  \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUDA driver error in line " << __LINE__;                   \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Error: ";                                                \
            errorStr << scai::common::cudaDriverErrorString( res );                 \
            errorStr << ", CUresult = " << res << "\n";                             \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw exception( errorStr.str() );                                      \
        }                                                                           \
    }

#define SCAI_CUDA_DRV_CALL_NOTHROW( call, msg )                                     \
    {                                                                               \
        CUresult res = call;                                                        \
        if ( CUDA_SUCCESS != res )                                                  \
        {                                                                           \
            std::cerr << "CUDA driver error in line " << __LINE__;                  \
            std::cerr << " of file " << __FILE__ << std::endl;                      \
            std::cerr << "  Msg  : " << msg << std::endl;                           \
            std::cerr << "  Call : " #call;                                         \
            std::cerr << "  Error: ";                                               \
            std::cerr << scai::common::cudaDriverErrorString( res );                \
            std::cerr << ", CUresult = " << res << "\n";                            \
            scai::common::Exception::addCallStack( std::cerr );                     \
        }                                                                           \
    }

#define SCAI_CUDA_DRV_CALL( call, msg )                                             \
    SCAI_CUDA_DRV_CALL_EXCEPTION( call, msg, scai::common::CUDAException )

#define SCAI_CUDA_RT_CALL(call, msg)                                                \
    {                                                                               \
        cudaError_t res = call;                                                     \
        if ( cudaSuccess != res )                                                   \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUDA runtime error in line " << __LINE__;                  \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Error: ";                                                \
            errorStr << cudaGetErrorString( res );                                  \
            errorStr << ", cudaError_t = " << res << "\n";                          \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw scai::common::Exception( errorStr.str() );                        \
        }                                                                           \
    }

#define SCAI_CUBLAS_CALL( call, msg )                                               \
    {                                                                               \
        cublasStatus_t res = call;                                                  \
        if ( CUBLAS_STATUS_SUCCESS != res )                                         \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUBLAS error in line " << __LINE__;                        \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Error: ";                                                \
            errorStr << scai::common::cublasErrorString( res );                     \
            errorStr << ", cublasStatus = " << res << "\n";                         \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw scai::common::CUDAException( errorStr.str() );                    \
        }                                                                           \
    }

#define SCAI_CUSPARSE_CALL( call, msg )                                             \
    {                                                                               \
        cusparseStatus_t res = call;                                                \
        if ( CUSPARSE_STATUS_SUCCESS != res )                                       \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUSparse error in line " << __LINE__;                      \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Error: ";                                                \
            errorStr << scai::common::cusparseErrorString( res );                   \
            errorStr << ", cusparseStatus = " << res << "\n";                       \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw scai::common::CUDAException( errorStr.str() );                    \
        }                                                                           \
    }

#define SCAI_CUFFT_CALL( call, msg )                                                \
    {                                                                               \
        cufftResult res = call;                                                     \
        if ( CUFFT_SUCCESS != res )                                                 \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUFFT error in line " << __LINE__;                         \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Error: ";                                                \
            errorStr << scai::common::cufftErrorString( res );                      \
            errorStr << ", cufftStatus = " << res << "\n";                          \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw scai::common::CUDAException( errorStr.str() );                    \
        }                                                                           \
    }

#if ( CUDART_VERSION >= 7050 )
#define SCAI_CUSOLVER_CALL( call, msg )                                             \
    {                                                                               \
        cusolverStatus_t res = call;                                                \
        if ( CUSOLVER_STATUS_SUCCESS != res )                                       \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "CUSolver error in line " << __LINE__;                      \
            errorStr << " of file " << __FILE__ << std::endl;                       \
            errorStr << "  Call : " #call;                                          \
            errorStr << "  Msg  : " << msg << std::endl;                            \
            errorStr << "  Error: ";                                                \
            errorStr << scai::common::cusolverErrorString( res );                   \
            errorStr << ", cusolverStatus = " << res << "\n";                       \
            scai::common::Exception::addCallStack( errorStr );                      \
            throw scai::common::CUDAException( errorStr.str() );                    \
        }                                                                           \
    }
#endif

#define SCAI_CHECK_CUDA_ACCESS                                                      \
    {                                                                               \
        CUcontext pctx;                                                             \
        SCAI_CUDA_DRV_CALL( cuCtxGetCurrent( &pctx ), "" );                         \
        SCAI_ASSERT( pctx, "No current context, forgotten SCAI_CONTEXT_ACCESS ?" )  \
    }

#define SCAI_CHECK_CUDA_ERROR                                                       \
    {                                                                               \
        SCAI_CUDA_RT_CALL( cudaGetLastError(), "last CUDA error" )                  \
    }

//    #define SCAI_CHECK_CUDA_ERROR


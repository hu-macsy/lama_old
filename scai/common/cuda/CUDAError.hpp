/**
 * @file CUDAError.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

// CUDA
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>

namespace scai
{

namespace common
{
/** Function that translates enum CUresult to strings. */

const char* cudaDriverErrorString( CUresult res );

/** Function that translates enum cublasStatus to strings. */

const char* cublasErrorString( cublasStatus_t res );

/** Function that translates enum cuparseStatus to strings. */

const char* cusparseErrorString( cusparseStatus_t res );

} /* end namespace common */

} /* end namespace scai */

/** Macro for CUDA driver API calls to catch errors */

#define SCAI_CUDA_DRV_CALL_EXCEPTION( call, msg, ExceptionClass )                       \
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
            errorStr << scai::common::cudaDriverErrorString( res );                     \
            errorStr << ", CUresult = " << res << "\n";                                 \
            scai::common::Exception::addCallStack( errorStr );                          \
            throw ExceptionClass( errorStr.str() );                                     \
        }                                                                               \
    }

#define SCAI_CUDA_DRV_CALL( call, msg )                                                 \
    SCAI_CUDA_DRV_CALL_EXCEPTION( call, msg, scai::common::Exception )                  

#define SCAI_CUDA_RT_CALL(call, msg)                                                    \
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
            scai::common::Exception::addCallStack( errorStr );                          \
            throw scai::common::Exception( errorStr.str() );                            \
        }                                                                               \
    }

#define SCAI_CUBLAS_CALL( call, msg )                                                   \
    {                                                                                   \
        cublasStatus_t res = call;                                                      \
        if ( CUBLAS_STATUS_SUCCESS != res )                                             \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "CUBLAS error in line " << __LINE__;                            \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Error: ";                                                    \
            errorStr << scai::common::cublasErrorString( res );                         \
            errorStr << ", cublasStatus = " << res << "\n";                             \
            scai::common::Exception::addCallStack( errorStr );                          \
            throw scai::common::Exception( errorStr.str() );                            \
        }                                                                               \
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
            throw scai::common::Exception( errorStr.str() );                        \
        }                                                                           \
    }

#define SCAI_CHECK_CUDA_ACCESS                                                          \
    {                                                                                   \
        CUcontext pctx;                                                                 \
        SCAI_CUDA_DRV_CALL( cuCtxGetCurrent( &pctx ), "" );                             \
        SCAI_ASSERT( pctx, "No current context, forgotten SCAI_CONTEXT_ACCESS ?" )    \
    }

#define SCAI_CHECK_CUDA_ERROR                                                         \
    {                                                                                 \
        SCAI_CUDA_RT_CALL( cudaGetLastError(), "last CUDA error" )                    \
    }

//    #define SCAI_CHECK_CUDA_ERROR


/**
 * @file CUDAAccess.cpp
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
 * @brief Implemenation of methods for class CUDAAccess.
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDACtx.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

namespace scai
{

namespace common
{

// The current CUDA device can be accessed globally, but should be thread-private

thread_local const CUDACtx* currentCUDACtx = NULL;

const CUDACtx* CUDAAccess::enable( const CUDACtx& ctx )
{
    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( ctx.getCUcontext() ), "could not push context" )
    const CUDACtx* last = currentCUDACtx;
    currentCUDACtx = &ctx;  // make it available globally in thread-private variable
    return last;
}

void CUDAAccess::disable( const CUDACtx* last )
{
    CUcontext tmp; // result variable for current context, not needed here
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
    currentCUDACtx = last;
    // last != NULL -> current context is last->getCUcontext()
    // last == NULL -> current context is 0
}

CUDAAccess::CUDAAccess( const CUDACtx& ctx ) : mCUcontext( ctx.getCUcontext() )
{
    mSaveDevice = enable( ctx );
}

CUDAAccess::~CUDAAccess()
{
    disable( mSaveDevice );
}

const CUDACtx& CUDAAccess::getCurrentCUDACtx()
{
    SCAI_ASSERT( currentCUDACtx, "Currently, no context is set" )
    return *currentCUDACtx;
}

} /* end namespace common */

} /* end namespace scai */

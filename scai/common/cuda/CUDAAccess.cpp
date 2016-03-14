/**
 * @file CUDAAccess.cpp
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
 * @brief Implemenation of methods for class CUDAAccess.
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDACtx.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Thread.hpp>

#include <iostream>

namespace scai
{

namespace common
{

// The current CUDA device can be accessed globally, but should be thread-private
// we can rely on the fact that thread-private variable is initialized with NULL 

static common::ThreadPrivatePtr<const CUDACtx> currentCUDACtx;

void CUDAAccess::enable( const CUDACtx& dev )
{
    // if ( currentCUDACtx.get() != NULL )
    // {
    //     std::cout << "WARN: enable of CUDA context while other context is set" << std::endl;
    // }

    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( dev.getCUcontext() ), "could not push context" )

    currentCUDACtx.set( &dev );
}

void CUDAAccess::disable()
{
    CUcontext tmp; // temporary for last context, not necessary to save it

    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )

    currentCUDACtx.set( NULL );
}

CUDAAccess::CUDAAccess( const CUDACtx& dev ) : mCUcontext( dev.getCUcontext() )
{
    mSaveDevice = currentCUDACtx.get();
    enable( dev );
}

CUDAAccess::~CUDAAccess()
{ 
    disable();

    currentCUDACtx.set( mSaveDevice );
}

const CUDACtx& CUDAAccess::getCurrentCUDACtx()
{
    const CUDACtx* current = currentCUDACtx.get();

    SCAI_ASSERT( current, "Currently, no device is accessed" )

    return *current;
}

} /* end namespace common */

} /* end namespace scai */

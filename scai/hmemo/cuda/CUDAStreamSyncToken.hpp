/**
 * @file CUDAStreamSyncToken.hpp
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
 * @brief Defintion of a SyncToken class that synchronizes with computations and
 *        memory transfers on CUDA devices.
 * @author Jiri Kraus, Thomas Brandes
 * @date 28.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/common/shared_ptr.hpp>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h> /* no diagnostic for this one */

// std
#include <list>
#include <memory>

namespace scai
{

namespace hmemo
{
    class CUDAContext;

    typedef common::shared_ptr<const CUDAContext> CUDAContextPtr;
}

namespace tasking
{

/** Class that sycnchronizes with a CUDA stream. */

class COMMON_DLL_IMPORTEXPORT CUDAStreamSyncToken: public SyncToken

{
public:

    /** Constructor for a sychronization token.
     *
     *  @param[in]  context  is the CUDAcontext of the stream
     *  @param[in]  stream   is the handle of the CUDA stream.
     *
     *  A pointer to the CUDA context is required to enable/disable it.
     */

    CUDAStreamSyncToken( hmemo::CUDAContextPtr context, CUstream stream );

    CUDAStreamSyncToken( hmemo::CUDAContextPtr context, CUstream stream, CUevent event );

    void setEvent( CUevent event )
    {
        mEvent = event;
    }

    virtual ~CUDAStreamSyncToken();

    virtual void wait();

    virtual bool probe() const;

    cudaStream_t getCUDAStream() const;

    void createTimingEvent( CUevent& event ) const;

    void createEvent( CUevent& event ) const;

    void getTime( float* time, CUevent& startEvent, CUevent& stopEvent ) const;

    bool probeEvent( const CUevent& stopEvent ) const;

    void recordEvent( const CUevent event );

    bool queryEvent( const CUevent event ) const;

    void synchronizeEvent( const CUevent event ) const;

private:

    hmemo::CUDAContextPtr mCUDAContext; // needed for synchronization

    const CUstream mStream;

    CUevent mEvent;
};

} /* end namespace hmemo */

} /* end namespace scai */

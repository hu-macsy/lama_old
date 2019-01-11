/**
 * @file CUDAStreamSyncToken.hpp
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
 * @brief Defintion of a SyncToken class that synchronizes with computations and
 *        memory transfers on CUDA devices.
 * @author Jiri Kraus, Thomas Brandes
 * @date 28.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// used classes

#include <scai/tasking/SyncToken.hpp>
#include <scai/tasking/cuda/CUDAStreamPool.hpp>
#include <scai/common/cuda/CUDACtx.hpp>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h> /* no diagnostic for this one */

// std
#include <list>
#include <memory>

namespace scai
{

namespace tasking
{

/** Class that sycnchronizes with a CUDA stream. */

class COMMON_DLL_IMPORTEXPORT CUDAStreamSyncToken: public SyncToken

{
public:

    /** Constructor for a sychronization token.
     *
     *  @param[in]  cuda   is the CUDAcontext of the stream
     *  @param[in]  type   is the type of the CUDA stream.
     *
     *  A pointer to the CUDA context is required to enable/disable it.
     */

    CUDAStreamSyncToken( const common::CUDACtx& cuda, const StreamType type );

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

    /** Get sync token in case of asynchronous execution should be started. */

    static CUDAStreamSyncToken* getCurrentSyncToken();

    /** Override SyncToken::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

private:

    const common::CUDACtx& mCUDA;   // needed for synchronization

    CUstream mStream;

    CUevent mEvent;
};

} /* end namespace hmemo */

} /* end namespace scai */

/**
 * @file CUDAStreamPool.hpp
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
 * @brief Pool to manage multiple streams for a CUDA context
 * @author Thomas Brandes
 * @date 08.03.2016
 */
#pragma once

#include <scai/common/config.hpp>

#include <scai/common/cuda/CUDACtx.hpp>

#include <scai/logging.hpp>

namespace scai
{

namespace tasking
{

enum class StreamType
{
    ComputeStream,
    TransferStream
};

/** For each CUDA device/context there will be some streams available
 *  that can be used for asynchronous computations and memory transfers.
 *
 *  The advantage of the pool is that create/destroy of the stream is only called once.
 */

class COMMON_DLL_IMPORTEXPORT CUDAStreamPool 

{
public:

    /** Get a stream of the pool, either for compute or memory transfer */

    CUstream reserveStream( StreamType type );

    /** Release a stream, no more used. */

    void releaseStream( CUstream stream );

    /** Get the stream pool for a CUDA context.
     *
     *  A new pool will might be created if not available.
     */

    static CUDAStreamPool& getPool( const common::CUDACtx& cuda );

    /** Release the stream pool for a CUDA context.
     *
     *  @throws an exception if not all streams have been released.
     */

    static void freePool( const common::CUDACtx& cuda );

    /** 
     *  Query if all streams are released 
     */
    bool isEmpty();

private:

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    const common::CUDACtx& mCUDA;

    /** Construct a pool of streams for a given CUDA device. */

    CUDAStreamPool( const common::CUDACtx& cuda );
    ~CUDAStreamPool();

    CUstream mTransferStream;
    CUstream mComputeStream;

    int mTransferReservations;
    int mComputeReservations;
};

} /* end namespace tasking */

} /* end namespace scai */

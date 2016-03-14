/**
 * @file CUDAStreamPool.hpp
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

struct streamtype
{
    typedef enum
    {
        ComputeStream,
        TransferStream

    } StreamType;
};

/** For each CUDA device/context there will be some streams available 
 *  that can be used for asynchronous computations and memory transfers.
 *
 *  The advantage of the pool is that create/destroy of the stream is only called once.
 */

class COMMON_DLL_IMPORTEXPORT CUDAStreamPool : public streamtype

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

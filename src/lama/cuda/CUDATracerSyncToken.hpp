/**
 * @file CUDATracerSyncToken.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief CUDATracerSyncToken.hpp
 * @author schubert
 * @date 09.11.2011
 * $Id$
 */
#ifndef LAMA_CUDATRACERSYNCTOKEN_HPP_
#define LAMA_CUDATRACERSYNCTOKEN_HPP_

#include <lama/SyncToken.hpp>

#include <lama/tracing/LAMABaseTracer.hpp>

#include <lama/cuda/CUDAStreamSyncToken.hpp>

#include <cuda.h>

/**
 * @brief lama::CUDATracerSyncToken
 */
class CUDATracerSyncToken: public lama::SyncToken
{
public:
    CUDATracerSyncToken( std::auto_ptr<LAMABaseTracer> tracer, lama::CUDAStreamSyncToken& cudaStreamSyncToken );

    virtual ~CUDATracerSyncToken();

    virtual void wait();

    virtual bool probe() const;

    void recordStopEvent();

private:

    CUDATracerSyncToken();

    CUDATracerSyncToken( const CUDATracerSyncToken& other );

    CUDATracerSyncToken& operator=( const CUDATracerSyncToken& other );

    CUevent mStartEvent;
    CUevent mStopEvent;

    lama::CUDAStreamSyncToken& mStreamSyncToken;

    std::auto_ptr<LAMABaseTracer> mTracer;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

#endif // LAMA_CUDATRACERSYNCTOKEN_HPP_

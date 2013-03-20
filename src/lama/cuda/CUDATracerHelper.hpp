/**
 * @file CUDATracerHelper.hpp
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
 * @brief CUDATracerHelper.hpp
 * @author schubert
 * @date 09.11.2011
 * $Id$
 */
#ifndef LAMA_CUDATRACERHELPER_HPP_
#define LAMA_CUDATRACERHELPER_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/tracing/LAMABaseTracer.hpp>
#include <lama/tracing/CUDATracerSyncToken.hpp>

#include <lama/cuda/CUDAStreamSyncToken.hpp>

/**
 * @brief CUDATracerHelper
 */
template<typename Tracer>
class LAMA_DLL_IMPORTEXPORT CUDATracerHelper
{
public:
    CUDATracerHelper( const char* name, const char* file, int lno, lama::CUDAStreamSyncToken& cudaStreamSyncToken );
    virtual ~CUDATracerHelper();

private:
    CUDATracerSyncToken* mTracerSyncToken;
};

template<typename Tracer>
CUDATracerHelper<Tracer>::CUDATracerHelper(
    const char* name,
    const char* file,
    int lno,
    lama::CUDAStreamSyncToken& cudaStreamSyncToken )
    : mTracerSyncToken( 0 )
{
    std::auto_ptr<LAMABaseTracer> tracer( new Tracer( name, file, lno ) );

    std::auto_ptr<CUDATracerSyncToken> cudaTracerSyncToken( new CUDATracerSyncToken( tracer, cudaStreamSyncToken ) );
    mTracerSyncToken = cudaTracerSyncToken.get();
    std::auto_ptr<lama::SyncToken> tracerSyncToken( cudaTracerSyncToken );
    cudaStreamSyncToken.pushSyncToken( tracerSyncToken );
}

template<typename Tracer>
CUDATracerHelper<Tracer>::~CUDATracerHelper()
{
    mTracerSyncToken->recordStopEvent();
}

#endif // LAMA_CUDATRACERHELPER_HPP_

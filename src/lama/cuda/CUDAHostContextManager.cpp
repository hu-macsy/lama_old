/**
 * @file CUDAHostContextManager.cpp
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
 * @brief Implementation of CUDA host manager class methods.
 * @author Thomas Brandes
 * @date 16.07.2011
 * $Id$
 */

// hpp
#include <lama/cuda/CUDAHostContextManager.hpp>

// others
#include <lama/cuda/CUDAHostContext.hpp>

#include <lama/exception/Exception.hpp>

namespace lama
{

// Explicit definition of static variable required.

CUDAHostContextManager CUDAHostContextManager::theInstance;

boost::weak_ptr<const CUDAHostContext> CUDAHostContextManager::contextInstance;

int CUDAHostContextManager::cudaDeviceNr = LAMA_DEFAULT_DEVICE_NUMBER;

/* ------------------------------------------------------------------------ */

CUDAHostContextManager::CUDAHostContextManager()
    : ContextManager( Context::Host )
{
    // This manager will NOT register itself as default, must be done explicitly
}

/* ------------------------------------------------------------------------ */

CUDAHostContextManager::~CUDAHostContextManager()
{
}

/* ------------------------------------------------------------------------ */

void CUDAHostContextManager::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "CUDAHostContextManager";
}

/* ------------------------------------------------------------------------ */

ContextPtr CUDAHostContextManager::getContext( int deviceNr )
{
    boost::shared_ptr<const CUDAHostContext> context;

    if ( deviceNr != LAMA_DEFAULT_DEVICE_NUMBER )
    {
        LAMA_LOG_WARN( logger, "Context number ignored for HostContext, deviceNr = " << deviceNr )
    }

    // use the last contextInstance if it is still valid

    if ( contextInstance.expired() )
    {
        // create a new instance of CUDAHostContext and keep it for further uses

        boost::shared_ptr<const CUDAContext> cudaContext = boost::dynamic_pointer_cast<const CUDAContext>(
                    ContextFactory::getContext( Context::CUDA, cudaDeviceNr ) );

        context = boost::shared_ptr<const CUDAHostContext>( new CUDAHostContext( cudaContext ) );

        contextInstance = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstance.lock();
    }

    return context;
}

void CUDAHostContextManager::setAsCurrent( ContextPtr cudaContext )
{
    if ( !cudaContext )
    {
        LAMA_THROWEXCEPTION( "no valid CUDAContext, is NULL" )
    }

    if ( cudaContext->getType() != Context::CUDA )
    {
        LAMA_THROWEXCEPTION( "no valid CUDAContext, is not CUDA" )
    }

    const CUDAContext& cc = static_cast<const CUDAContext&>( *cudaContext );

    // we will just store the device number so that we can easily get it back

    cudaDeviceNr = cc.getDeviceNr();

    // set singleton object for CUDAHost context manager in ContextFactory

    theInstance.registerFactory();
}

} // namespace lama


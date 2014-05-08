/**
 * @file CUDAawareMPICommunicatorManager.cpp
 *
 * @license
 * Copyright (c) 2014
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
 * @brief CUDAawareMPICommunicatorManager.cpp
 * @author lschubert
 * @date 11.03.2014
 * @since 1.1.0
 */

// hpp
#include <lama/mpi/CUDAawareMPICommunicatorManager.hpp>

// others
#include <lama/CommunicatorFactory.hpp>

#include <cstdio>

using namespace std;

#define COMMUNICATOR_TYPE "CUDAawareMPI"

namespace lama
{

/* ----  Static variables  -------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( CUDAawareMPICommunicatorManager::logger, "CommunicatorManager.CUDAawareMPICommunicatorManager" )

// make sure that static initialization is called.

bool CUDAawareMPICommunicatorManager::__init = CUDAawareMPICommunicatorManager::init();

// static initialization adds an MPI communicator manager to the communicator factory

bool CUDAawareMPICommunicatorManager::init()
{
    // logger should already be available as it is initialized before __init

    LAMA_LOG_DEBUG( logger, "create MPI communicator manager and add it to CommunicatorFactory" )

    boost::shared_ptr<CommunicatorManager> manager( new CUDAawareMPICommunicatorManager() );
    CommunicatorFactory::getFactory().addCommunicatorManager( COMMUNICATOR_TYPE, manager );
    return true;
}

/* -------------------------------------------------------------------------- */

CUDAawareMPICommunicatorManager::CUDAawareMPICommunicatorManager()
    : CommunicatorManager( COMMUNICATOR_TYPE )
{
}

/* -------------------------------------------------------------------------- */

CUDAawareMPICommunicatorManager::~CUDAawareMPICommunicatorManager()
{
    LAMA_LOG_INFO( logger, "~CUDAawareMPICommunicatorManager(), use count MPI communicator = " << mCommInstance.use_count() )

    if ( mCommInstance.use_count() > 1 )
    {
        LAMA_LOG_WARN( logger,
                       "MPICommunicator has " << mCommInstance.use_count() - 1 << " remaining references, seems that not all LAMA data structures have been freed" )
    }

    // Note: destructor will free commInstance and therefore call MPI_Finalize
}

/* -------------------------------------------------------------------------- */

CommunicatorPtr CUDAawareMPICommunicatorManager::getCommunicator( int& argc, char**& argv )
{
    if ( !mCommInstance )
    {
        LAMA_LOG_INFO( logger, "create new MPICommunicator" )

        // create a new instance of MPICommunicator, will call MPI_Init

        mCommInstance = boost::shared_ptr<const CUDAawareMPICommunicator>( new CUDAawareMPICommunicator( argc, argv ) );
    }

    return mCommInstance;
}

}

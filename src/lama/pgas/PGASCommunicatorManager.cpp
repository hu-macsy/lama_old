/**
 * @file PGASCommunicatorManager.cpp
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
 * @brief PGASCommunicatorManager.cpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASCommunicatorManager.hpp>
#include <lama/pgas/PGASCommunicator.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <boost/shared_ptr.hpp>

#define COMMUNICATOR_TYPE "PGAS"

namespace lama
{

/* ----  Static variables  -------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( PGASCommunicatorManager::logger, "CommunicatorManager.PGASCommunicatorManager" );

// make sure that static initialization is called.

bool PGASCommunicatorManager::__init = PGASCommunicatorManager::init();

// static initialization adds an PGAS communicator manager to the communicator factory

bool PGASCommunicatorManager::init()
{
    // logger should already be available as it is initialized before __init
    LAMA_LOG_DEBUG( logger, "create PGAS communicator manager and add it to CommunicatorFactory" );
    boost::shared_ptr<CommunicatorManager> manager( new PGASCommunicatorManager() );
    CommunicatorFactory::getFactory().addCommunicatorManager( COMMUNICATOR_TYPE, manager );
    PGASInterface::getInstance(); //call get Instance so the PGASInterface gets initialized
    return true;
}

/* -------------------------------------------------------------------------- */

PGASCommunicatorManager::PGASCommunicatorManager()
    : CommunicatorManager( COMMUNICATOR_TYPE )
{
}

/* -------------------------------------------------------------------------- */

PGASCommunicatorManager::~PGASCommunicatorManager()
{
    LAMA_LOG_INFO( logger, "~PGASCommunicatorManager(), use count PGAS communicator = " << mCommInstance.use_count() );

    if( mCommInstance.use_count() > 1 )
    {
        LAMA_LOG_WARN( logger,
                       "PGASCommunicator has " << mCommInstance.use_count() - 1 << " remaining references, seems that not all LAMA data structures have been freed" );
    }
}

/* -------------------------------------------------------------------------- */

CommunicatorPtr PGASCommunicatorManager::getCommunicator( int& argc, char**& argv )
{
    if( !mCommInstance )
    {
        LAMA_LOG_INFO( logger, "create new PGASCommunicator" );

        // create a new instance of PGASCommunicator, will call PGAS_Init

        mCommInstance = boost::shared_ptr<const PGASCommunicator>( new PGASCommunicator( argc, argv ) );
    }

    return mCommInstance;
}

}

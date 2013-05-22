/**
 * @file CyclicDistributionManager.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implemenation of class lama::CyclicDistributionManager
 * @author Thomas Brandes
 * @date 20.12.2012
 * $Id$
 */

// hpp
#include <lama/distribution/CyclicDistributionManager.hpp>

// others
#include <lama/distribution/CyclicDistribution.hpp>

#include <lama/CommunicatorFactory.hpp>
#include <lama/DistributionFactory.hpp>

#include <cstdio>

using namespace std;

#define DISTRIBUTION_TYPE "CYCLIC"

namespace lama
{

/* ----  Static variables  -------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( CyclicDistributionManager::logger, "DistributionManager.CyclicDistributionManager" )

// make sure that static initialization is called.

bool CyclicDistributionManager::__init = CyclicDistributionManager::init();

// static initialization adds an NO communicator manager to the communicator factory

bool CyclicDistributionManager::init()
{
    // logger should already be available as it is initialized before __init

    LAMA_LOG_DEBUG( logger, "create NO communicator manager and add it to DistributionFactory" )

    boost::shared_ptr<DistributionManager> manager( new CyclicDistributionManager() );
    DistributionFactory::getFactory().addDistributionManager( DISTRIBUTION_TYPE, manager );
    return true;
}

/* -------------------------------------------------------------------------- */

CyclicDistributionManager::CyclicDistributionManager()
    : DistributionManager( DISTRIBUTION_TYPE )
{
}

/* -------------------------------------------------------------------------- */

CyclicDistributionManager::~CyclicDistributionManager()
{
    LAMA_LOG_DEBUG( logger, "~CyclicDistributionManager" )
}

/* -------------------------------------------------------------------------- */

DistributionPtr CyclicDistributionManager::getDistribution( const IndexType size, const std::vector<int>& arguments )
{
    int block_size = size; // default

    if ( arguments.size() == 0 )
    {
        LAMA_LOG_WARN( logger, "argument list should have exactly one int value" )
    }
    else
    {
        block_size = arguments[0];
    }

    CommunicatorPtr comm = CommunicatorFactory::get(); // default communicator

    return boost::shared_ptr<const CyclicDistribution>( new CyclicDistribution( size, block_size, comm ) );
}

}

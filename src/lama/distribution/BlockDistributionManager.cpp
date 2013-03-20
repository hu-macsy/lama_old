/**
 * @file BlockDistributionManager.cpp
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
 * @brief Implemenation of class lama::BlockDistributionManager
 * @author Thomas Brandes
 * @date 20.12.2012
 * $Id$
 */

// hpp
#include <lama/distribution/BlockDistributionManager.hpp>
#include <lama/CommunicatorFactory.hpp>

// others
#include <lama/DistributionFactory.hpp>

#include <cstdio>

using namespace std;

#define DISTRIBUTION_TYPE "BLOCK"

namespace lama
{

/* ----  Static variables  -------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( BlockDistributionManager::logger, "DistributionManager.BlockDistributionManager" );

// make sure that static initialization is called.

bool BlockDistributionManager::__init = BlockDistributionManager::init();

// static initialization adds an NO communicator manager to the communicator factory

bool BlockDistributionManager::init()
{
    // logger should already be available as it is initialized before __init

    LAMA_LOG_DEBUG( logger, "create NO communicator manager and add it to DistributionFactory" );

    boost::shared_ptr<DistributionManager> manager( new BlockDistributionManager() );
    DistributionFactory::getFactory().addDistributionManager( DISTRIBUTION_TYPE, manager );
    return true;
}

/* -------------------------------------------------------------------------- */

BlockDistributionManager::BlockDistributionManager()
    : DistributionManager( DISTRIBUTION_TYPE )
{
}

/* -------------------------------------------------------------------------- */

BlockDistributionManager::~BlockDistributionManager()
{
    LAMA_LOG_INFO( logger, "~BlockDistributionManager" );
}

/* -------------------------------------------------------------------------- */

DistributionPtr BlockDistributionManager::getDistribution( const IndexType size, const std::vector<int>& )
{
    // ToDo: check for no additional arguments

    CommunicatorPtr comm = CommunicatorFactory::get(); // default communicator

    return boost::shared_ptr<const BlockDistribution>( new BlockDistribution( size, comm ) );
}

}

/**
 * @file NoDistributionManager.cpp
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
 * @brief Implemenation of class lama::NoDistributionManager
 * @author Thomas Brandes
 * @date 20.12.2012
 * @since 1.0.0
 */

// hpp
#include <lama/distribution/NoDistributionManager.hpp>

// others
#include <lama/distribution/NoDistribution.hpp>

#include <lama/DistributionFactory.hpp>

#include <cstdio>

using namespace std;

#define DISTRIBUTION_TYPE "NO"

namespace lama
{

/* ----  Static variables  -------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( NoDistributionManager::logger, "DistributionManager.NoDistributionManager" )

// make sure that static initialization is called.

bool NoDistributionManager::__init = NoDistributionManager::init();

// static initialization adds an NO communicator manager to the communicator factory

bool NoDistributionManager::init()
{
    // logger should already be available as it is initialized before __init

    LAMA_LOG_DEBUG( logger, "create NO distribution manager and add it to DistributionFactory" )

    boost::shared_ptr<DistributionManager> manager( new NoDistributionManager() );
    DistributionFactory::getFactory().addDistributionManager( DISTRIBUTION_TYPE, manager );
    return true;
}

/* -------------------------------------------------------------------------- */

NoDistributionManager::NoDistributionManager()
    : DistributionManager( DISTRIBUTION_TYPE )
{
}

/* -------------------------------------------------------------------------- */

NoDistributionManager::~NoDistributionManager()
{
    LAMA_LOG_INFO( logger, "~NoDistributionManager" )
}

/* -------------------------------------------------------------------------- */

DistributionPtr NoDistributionManager::getDistribution( const IndexType size, const std::vector<int>& arguments )
{
    if ( arguments.size() > 0 )
    {
        LAMA_LOG_WARN( logger, "additional arguments are ignored" )
    }

    return boost::shared_ptr<const NoDistribution>( new NoDistribution( size ) );
}

}

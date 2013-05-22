/**
 * @file CyclicDistributionManager.hpp
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
 * @brief Definition of manager class for block distributions.
 * @author Thomas Brandes
 * @date 20.10.2012
 * $Id$
 */

#ifndef LAMA_CYCLIC_DISTRIBUTION_MANAGER_HPP_
#define LAMA_CYCLIC_DISTRIBUTION_MANAGER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/DistributionManager.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** Singleton class to manage replicated distributions (instance on each processor).
 */

class CyclicDistributionManager: public DistributionManager
{

public:

    /** Provide implementation for the abstract routine of class DistributionManager. */

    virtual DistributionPtr getDistribution( const IndexType size, const std::vector<int>& arguments );

    virtual ~CyclicDistributionManager();

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static bool __init; //!< used to force call of init()

    static bool init();//!< create an instance that is added to distribution factory

    CyclicDistributionManager();
};

}

#endif // LAMA_CYCLIC_DISTRIBUTION_MANAGER_HPP_

/**
 * @file DistributionManager.hpp
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
 * @brief Base class for a distribution manager that is used for the distribution factory.
 * @author Thomas Brandes
 * @date 23.02.2011
 * $Id$
 */

#ifndef LAMA_DISTRIBUTION_MANAGER_HPP_
#define LAMA_DISTRIBUTION_MANAGER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/distribution/Distribution.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** @brief Base class for a distribution manager.
 *
 *  For each distribution type a distribution manager must be available.
 *  The distribution manager registers in the DistributionFactory and
 *  delivers a shared pointer for a distribution.
 *
 */

class DistributionManager: private NonCopyable
{
public:

    virtual ~DistributionManager();

    /** @brief Method that returns a distribution.
     *
     * @param[in] size is the number of elements to be distributed
     * @param[in] arguments is a vector with additional arguments
     *
     */
    virtual DistributionPtr getDistribution( const IndexType size, const std::vector<int>& arguments ) = 0;

protected:

    /** Constructor of a DistributionManager, type must be specified. */

    DistributionManager( const char* type );

    std::string mDistributionType; //!< type of distribution managed

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_DISTRIBUTION_MANAGER_HPP_

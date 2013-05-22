/**
 * @file CommunicatorManager.hpp
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
 * @brief CommunicatorManager.hpp
 * @author Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */

#ifndef LAMA_COMMUNICATOR_MANAGER_HPP_
#define LAMA_COMMUNICATOR_MANAGER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/Communicator.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** @brief Base class for a communicator manager.
 *
 *  For each communicator type a communicator manager must be available.
 *  The communicator manager registers in the CommunicatorFactory and
 *  delivers a shared pointer for a communicator.
 *
 */

class CommunicatorManager: private NonCopyable
{
public:

    virtual ~CommunicatorManager();

    /** @brief Method that returns a communicator.
     *
     * @param[in,out]  argc number of arguments
     * @param[in,out]  argv array of command line arguments
     *
     */
    virtual CommunicatorPtr getCommunicator( int& argc, char** & argv ) = 0;

protected:

    /** Constructor of a CommunicatorManager, type must be specified. */

    CommunicatorManager( const char* type );

    std::string mCommunicatorType; //!< type of communicator managed

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_COMMUNICATOR_MANAGER_HPP_

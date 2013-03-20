/**
 * @file PGASCommunicatorManager.hpp
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
 * @brief PGASCommunicatorManager.hpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#ifndef LAMA_PGAS_COMMUNICATOR_MANAGER_HPP_
#define LAMA_PGAS_COMMUNICATOR_MANAGER_HPP_

#include <lama/CommunicatorManager.hpp>
#include <lama/pgas/PGASCommunicator.hpp>

#include <logging/logging.hpp>

namespace lama
{

/** Singleton class to manage PGAS Communicator.
 *
 *  Regarding PGAS there is the restriction that after an PGAS_FINALIZE
 *  another call of PGAS_INIT is not allowed.
 */

class PGASCommunicatorManager: public CommunicatorManager
{

public:

    virtual CommunicatorPtr getCommunicator( int& argc, char** & argv );

    virtual ~PGASCommunicatorManager();

private:

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

    static bool __init; //!< used to force call of init()

    static bool init();//!< create an instance that is added to communicator factory

    PGASCommunicatorManager();

    /** PGAS communicator manager keeps a shared pointer to avoid multiple calls
     *  of PGAS_Init and PGAS_Finalize.
     */

    boost::shared_ptr<const PGASCommunicator> mCommInstance;
};

}

#endif // LAMA_PGAS_COMMUNICATOR_MANAGER_HPP_

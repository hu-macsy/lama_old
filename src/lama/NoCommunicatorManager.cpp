/**
 * @file NoCommunicatorManager.cpp
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
 * @brief Implementation of methods for manager class of NoCommunicator.
 * @author Thomas Brandes
 * @date 15.03.2011
 * @since 1.0.0
 */

// hpp
#include <lama/NoCommunicatorManager.hpp>

// others
#include <lama/NoCommunicator.hpp>
#include <lama/CommunicatorFactory.hpp>

// boost
#include <boost/shared_ptr.hpp>

using namespace std;

namespace lama
{

// make sure that static initialization is called.

bool NoCommunicatorManager::__init = NoCommunicatorManager::init();

bool NoCommunicatorManager::init()
{
    // static initialization adds a no communicator manager to the communicator factory

    boost::shared_ptr<CommunicatorManager> manager( new NoCommunicatorManager() );

    CommunicatorFactory::getFactory().addCommunicatorManager( "none", manager );

    return true;
}

NoCommunicatorManager::NoCommunicatorManager()
    : CommunicatorManager( "none" )
{
}

NoCommunicatorManager::~NoCommunicatorManager()
{
}

CommunicatorPtr NoCommunicatorManager::getCommunicator( int& /* argc */, char** & /* argv */)
{
    boost::shared_ptr<NoCommunicator> communicator;

    // use the last communicatorInstance if it is still valid

    if ( mCommunicatorInstance.expired() )
    {
        // create a new instance of NoCommunicator and keep it for further uses

        communicator = boost::shared_ptr<NoCommunicator>( new NoCommunicator() );

        mCommunicatorInstance = communicator;
    }
    else
    {
        // the last communicator instance is still valid, so we return new shared pointer to it

        communicator = mCommunicatorInstance.lock();
    }

    return communicator;
}

}


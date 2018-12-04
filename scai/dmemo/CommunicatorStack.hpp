/**
 * @file CommunicatorStack.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Class for stack of communicators used for tasking       
 * @author Thomas Brandes
 * @date 15.11.2018
 */
#pragma once

#include <scai/dmemo/Communicator.hpp>
#include <stack>

namespace scai
{

namespace dmemo
{

/** Stack used for task parallelism to use a sub-communicator in a scope. */

class CommunicatorStack 
{
public:

    static void push( CommunicatorPtr comm )
    {
        instance.push( comm );
    }

    static void pop()
    {
        instance.pop();
    }

    static bool empty()
    {
        return instance.empty();
    }

    static CommunicatorPtr top()
    {
        return instance.top();
    }

private:
 
    CommunicatorStack();

    static std::stack<CommunicatorPtr> instance;
};

class ScopedCommunicatorRecord 
{
public:

    ScopedCommunicatorRecord( CommunicatorPtr comm )
    {
        CommunicatorStack::push( comm );
    }

    ~ScopedCommunicatorRecord()
    {
        CommunicatorStack::pop();
    }
};

}

}

/* -------------------------------------------------------------------------- */

/** Macro that defines a new current communicator for the actual scope
 *
 *  \code
 *      auto comm = Communicator::getCommunicatorPtr();
 *      PartitionId color = ...;
 *      auto commTask = comm->split( color );
 *      {
 *          SCAI_DMEMO_TASK( commTask )
 *          auto dist = std::make_shared<BlockDistribution>( n );  // distributes onto commTask
 *          ...
 *      }
 *  \endcode
 */
#define SCAI_DMEMO_TASK( comm ) scai::dmemo::ScopedCommunicatorRecord SCAI_Comm_( comm );


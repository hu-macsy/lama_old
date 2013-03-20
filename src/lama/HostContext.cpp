/**
 * @file HostContext.cpp
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
 * @brief HostContext.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 * $Id$
 */

// hpp
#include <lama/HostContext.hpp>

// others
#include <lama/DefaultHostContext.hpp>

#include <lama/task/TaskSyncToken.hpp>

#include <cstdio>

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Constructors / Destructors                                               */
/* ------------------------------------------------------------------------- */

HostContext::HostContext()
    : Context( Host )
{
}

HostContext::~HostContext()
{
}

/* ------------------------------------------------------------------------- */

bool HostContext::canUseData( const Context& other ) const
{
    // same object by pointer can always use same data.

    if ( this == &other )
    {
        return true;
    }

    // different Host devices can use same data

    if ( other.getType() == Host )
    {
        return true;
        // equal if other is HostContext and has same host type
        // const HostContext& otherHost = static_cast<const HostContext&> (other);
        // return otherHost.getHostType() == getHostType();
    }

    return false;
}

/* ------------------------------------------------------------------------- */

void HostContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "HostContext";
}

std::auto_ptr<SyncToken> HostContext::getSyncToken() const
{
    // on Host we will run asynchronous computations as a task

    return std::auto_ptr<SyncToken>( new TaskSyncToken() );
}

}

/**
 * @file PGASSyncToken.cpp
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
 * @brief Implementation of methods for class PGASSyncToken.
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASSyncToken.hpp>
#include <lama/pgas/PGASInterface.hpp>

namespace lama
{

PGASSyncToken::PGASSyncToken()
{
    LAMA_LOG_INFO( logger, "PGASSyncToken constructed" )

}

PGASSyncToken::~PGASSyncToken()
{
    if ( !isSynchronized() )
    {
        LAMA_LOG_DEBUG( logger, *this << ": synchnronized at destructor" )
        wait();
    }
}

void PGASSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "PGASSyncToken, synchronized = " << isSynchronized() << " )";
}

void PGASSyncToken::wait()
{
    if ( isSynchronized() )
    {
        LAMA_LOG_WARN( logger, *this << ": waiting twice" )

        return; // do not wait twice, especially do not clean-up twice
    }

    LAMA_LOG_INFO( logger, *this << ": wait" )

    PGASInterface::getInstance()->syncronizeAll();

    LAMA_LOG_INFO( logger, *this << ": synchronized, clean up and free accesses" )

    setSynchronized();
}

bool PGASSyncToken::probe() const
{
    LAMA_LOG_WARN( logger, *this << ": probing is not possible --> probe has the same effect as isSyncronized()" )
    return isSynchronized();
}

}

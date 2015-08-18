/**
 * @file ContextAccess.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Implemenations for class ContextAccess.
 *
 * @author Thomas Brandes
 * @date 14.07.2011
 */

// hpp

#include <scai/hmemo/ContextAccess.hpp>

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( ContextAccess::logger, "ContextAccess" )

ContextAccess::ContextAccess( ContextPtr context, const char* file, int line )
                : mContext( *context ), mReleased( false ), mFile( file ), mLine( line )
{
    SCAI_LOG_INFO( logger, *this << " enabled" )

    mContext.enable( mFile, mLine );
}

void ContextAccess::release()
{
    if( mReleased )
    {
        return;
    }

    SCAI_LOG_INFO( logger, *this << " released" )

    mContext.disable( mFile, mLine );
    mReleased = false;
}

ContextAccess::~ContextAccess()
{
    release();
}

void ContextAccess::writeAt( std::ostream& stream ) const
{
    stream << "Access of " << mContext;
    stream << " at " << mFile << "( line = " << mLine << " )";
}

} /* end namespace hmemo */

} /* end namespace scai */

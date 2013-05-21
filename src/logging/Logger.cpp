/**
 * @file Logger.cpp
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
 * @brief Logger.cpp
 * @author Thomas Brandes
 * @date 01.03.2011
 * $Id$
 */

// hpp
#include <logging/Logger.hpp>

namespace log4lama
{

Logger::Logger( const std::string& name, Logger* const parent )
    : mName( name ), mSetFlag( false ), mLevel( WARN ), mParent( parent )
{
    if ( mParent != NULL )
    {
        mParent->mSons.push_back( this );
        mLevel = mParent->getLevel();
    }
}

Logger::~Logger()
{
}

std::string Logger::getFullName() const
{
    std::string fullname = "";

    if ( mParent == NULL )
    {
        fullname = "<root>";
    }
    else if ( mParent->isRootLogger() )
    {
        fullname = mName;
    }
    else
    {
        fullname = mParent->getFullName() + "." + mName;
    }

    return fullname;
}

Level Logger::getEffectiveLevel() const
{
    Level level = mLevel;

    if ( !mSetFlag && mParent != NULL )
    {
        level = mParent->getEffectiveLevel();
    }

    return level;
}

void Logger::setLevel( const Level level, const bool force/*= true*/)
{
    // if level is already set and there is no force => return
    if ( !force && mSetFlag )
    {
        return;
    }

    mLevel = level;

    mSetFlag = force;

    // traverse the sons but do no longer force
    for ( size_t i = 0; i < mSons.size(); i++ )
    {
        mSons[i]->setLevel( level, false );
    }
}

bool Logger::isRootLogger() const
{
    return mParent == NULL;
}

} // namespace log4lama

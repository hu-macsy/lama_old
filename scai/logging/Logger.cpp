/**
 * @file Logger.cpp
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
 * @brief Implementation of methods for Logger to deal with name hierarchies.
 * @author Thomas Brandes
 * @date 01.03.2011
 */

// hpp
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace logging
{

Logger::Logger( const std::string& name, Logger* const parent )
    : mName( name ), mSetFlag( false ), mLevel( level::WARN ), mParent( parent )
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

level::Level Logger::getEffectiveLevel() const
{
    level::Level level = mLevel;

    if ( !mSetFlag && mParent != NULL )
    {
        level = mParent->getEffectiveLevel();
    }

    return level;
}

void Logger::setLevel( const level::Level level, const bool force/*= true*/ )
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

} /* end namespace logging */

} /* end namespace scai */

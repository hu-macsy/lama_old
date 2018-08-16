/**
 * @file LoggerProvider.cpp
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
 * @brief Implementation of methods for class LoggerProvider.
 * @author Thomas Brandes
 * @date 02.03.2011
 */

// hpp
#include <scai/logging/LoggerProvider.hpp>

// local library
#include <scai/logging/Logger.hpp>

// std
#include <stdexcept>
#include <sstream>

namespace scai
{

namespace logging
{

// The LoggerProvier will be created at first access.

LoggerProvider* LoggerProvider::theProvider = NULL;

LoggerProvider& LoggerProvider::getProvider()
{
    if ( !theProvider )
    {
        theProvider = new LoggerProvider();
    }

    return *theProvider;
}

LoggerProvider::LoggerProvider()
    : mLoggerCreator( 0 )
{
}

LoggerProvider::~LoggerProvider()
{
}

void LoggerProvider::setLoggerCreator( AbstractLoggerCreator* creator )
{
    mLoggerCreator = creator;
}

Logger& LoggerProvider::getInstance( const std::string& name ) const
{
    if ( !mLoggerCreator )
    {
        mLoggerCreator = &theLoggerCreator();
    }

    std::vector<std::string> tokens;
    // Skip delimiters at beginning.
    std::string::size_type lastPos = name.find_first_not_of( ".", 0 );
    // Find first "non-delimiter".
    std::string::size_type pos = name.find_first_of( ".", lastPos );

    while ( std::string::npos != pos || std::string::npos != lastPos )
    {
        // Found a token, add it to the vector.
        tokens.push_back( name.substr( lastPos, pos - lastPos ) );
        // Skip delimiters.  Note the "not_of"
        lastPos = name.find_first_not_of( ".", pos );
        // Find next "non-delimiter"
        pos = name.find_first_of( ".", lastPos );
    }

    if ( !mLoggerCreator )
    {
        std::ostringstream errorMsg;
        errorMsg << "no LoggerCreator set in LoggerProvider" << ", cannot get instance for " << name;
        throw std::runtime_error( errorMsg.str() );
    }

    // now find the logger in the hierarchy tree
    Logger* instance = &mLoggerCreator->getRoot();

    for ( size_t i = 0; i < tokens.size(); i++ )
    {
        Logger* son = NULL; // will point to son that corresponds next token

        for ( size_t s = 0; s < instance->mSons.size(); s++ )
        {
            Logger* candidate = instance->mSons[s];

            if ( candidate->mName == tokens[i] )
            {
                son = candidate;
                break;
            }
        }

        // create a new son if not found

        if ( son == NULL )
        {
            // logger not available, so create it
            son = mLoggerCreator->create( tokens[i], instance );
        }

        // go to the next deeper level
        instance = son;
    }

    return *instance;
}

} /* end namespace logging */

} /* end namespace scai */

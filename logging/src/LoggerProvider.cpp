/**
 * @file LoggerProvider.cpp
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
 * @brief Implementation of methods for class LoggerProvider.
 * @author Thomas Brandes
 * @date 02.03.2011
 */

// hpp
#include <logging/LoggerProvider.hpp>

// others
#include <logging/Logger.hpp>

#include <stdexcept>

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

} // namespace logging

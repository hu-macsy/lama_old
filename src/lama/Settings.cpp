/**
 * @file Settings.cpp
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
 * @brief Managing some settings for CUDA specified by environment variables
 * @author Thomas Brandes
 * @date 04.05.2013
 * @since 1.0.0
 */

#include <lama/Settings.hpp>
#include <lama/Communicator.hpp>

#include <cstdio>
#include <cstring>
#include <cstdlib>

extern char **environ;

namespace lama
{

LAMA_LOG_DEF_LOGGER( Settings::logger, "Settings" )

/* ----------------------------------------------------------------------------- */

bool Settings::convertValue( int& flag, const char* stringVal )
{
    int nread = sscanf( stringVal, "%d", &flag );

    if( nread == 1 )
    {
        return true;
    }

    return false;
}

bool Settings::convertYesNoString( bool& flag, const char* stringVal )
{
    char key = static_cast<char>( toupper( stringVal[0] ) );

    bool done = true; // becomes false if no legal value has been found

    // to upper

    if( key == '0' )
    {
        flag = false;
    }
    else if( key == '1' )
    {
        flag = true;
    }
    else if( key == 'J' )
    {
        flag = true;
    }
    else if( key == 'Y' )
    {
        flag = true;
    }
    else if( key == 'T' )
    {
        flag = true;
    }
    else if( key == 'N' )
    {
        flag = false;
    }
    else if( key == 'F' )
    {
        flag = false;
    }
    else
    {
        // could not identify meaning
        done = false;
    }

    return done;
}

/* ----------------------------------------------------------------------------- */

bool Settings::getEnvironment( bool& flag, const char* envVarName )
{
    const char* env = getenv( envVarName );

    if( !env )
    {
        LAMA_LOG_INFO( logger, envVarName << " not set, will use other default" )

        return false; // no initialization by environment
    }

    bool done = convertYesNoString( flag, env );

    if( !done )
    {
        LAMA_LOG_ERROR( logger,
                        "Environment variable " << envVarName << "=" << env << ", is illegal setting, assume FALSE" )

        flag = false;
    }

    return true; // environment variable was available
}

bool Settings::getEnvironment( int& val, const char* envVarName )
{
    const char* env = getenv( envVarName );

    if( !env )
    {
        LAMA_LOG_INFO( logger, envVarName << " not set, will select by compute capability" )

        return false; // no initialization by environment
    }

    bool done = convertValue( val, env );

    if( !done )
    {
        LAMA_LOG_ERROR( logger,
                        "Environment variable " << envVarName << "=" << env << ", is illegal setting, assume FALSE" )

        return false;
    }

    return true; // environment variable was available
}

bool Settings::getEnvironment( std::string& val, const char* envVarName )
{
    const char *env = getenv( envVarName );

    if( env )
    {
        val = env;

        LAMA_LOG_INFO( logger, envVarName << " = " << val );

        return true;
    }

    LAMA_LOG_INFO( logger, envVarName << " not set" );

    return false;
}

bool Settings::getEnvironment( std::string& val, const char* envVarName, const Communicator& comm )
{
    bool isRoot = comm.getRank() == 0;

    bool hasSet = false;

    if( isRoot )
    {
        hasSet = getEnvironment( val, envVarName );

        if( hasSet )
        {
            comm.bcast( val, 0 );
        }
        else
        {
            std::string dummy = "";
            comm.bcast( dummy, 0 );
        }
    }
    else
    {
        std::string rootVal;

        comm.bcast( rootVal, 0 );

        bool hasRootSet = rootVal.length() > 0;

        if( hasRootSet )
        {
            val = rootVal;

            hasSet = true;

            // if process has own environment variable, overwrite it

            getEnvironment( val, envVarName );
        }
        else
        {
            hasSet = getEnvironment( val, envVarName );
        }
    }

    LAMA_LOG_INFO( logger, comm << ": " << envVarName << "=" << val );

    return hasSet;
}

bool Settings::init()
{
    int i = 0;
    char *s = *environ;

    for( ; s; i++ )
    {
        if( strncmp( s, "LAMA_", 5 ) == 0 )
        {
            printf( "%s\n", s );
        }

        s = *( environ + i );
    }

    return true;
}

} // namespace

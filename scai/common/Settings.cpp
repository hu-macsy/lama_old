/**
 * @file Settings.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Managing some settings for CUDA specified by environment variables
 * @author Thomas Brandes
 * @date 16.01.2016
 */

// hpp
#include <scai/common/Settings.hpp>

// std
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iostream>

extern char **environ;

namespace scai
{

namespace common
{

/* ----------------------------------------------------------------------------- */

void Settings::parseArgs( int& argc, const char* argv[] )
{
    const bool replace = true;   // command line args overwrite environment settings

    int unused_args = 0;

    for ( int i = 0; i < argc; ++i )
    {
        if ( strncmp( argv[i], "--SCAI_", 7 ) == 0 )
        {
            std::string arg = argv[i] + 2;

            std::string::size_type equalPos = arg.find_first_of( "=", 0 );

            if ( std::string::npos != equalPos )
            {
                // tokenize it  name = val

                std::string name = arg.substr( 0, equalPos );
                std::string val  = arg.substr( equalPos + 1 );

                putEnvironment( name.c_str(), val.c_str(), replace );
            }
        }
        else
        {
            argv[unused_args++] = argv[i];
        }
    }

    argc = unused_args;
}

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
    std::string val;

    if ( !getEnvironment( val, envVarName ) )
    {
        return false;  // environment variable not set
    }

    bool done = convertYesNoString( flag, val.c_str() );

    if ( !done )
    {
        return false;  //  no boolean value
    }

    return true; // environment variable was available
}

/* ----------------------------------------------------------------------------- */

bool Settings::getEnvironment( int& val, const char* envVarName )
{
    std::string envVal;

    if ( !getEnvironment( envVal, envVarName ) )
    {
        return false; // no initialization by environment
    }

    bool done = convertValue( val, envVal.c_str() );

    if ( !done )
    {
        return false;
    }

    return true; // environment variable was available
}

/* ----------------------------------------------------------------------------- */

bool Settings::getEnvironment( std::string& val, const char* envVarName )
{
    const char* env = getenv( envVarName );

    if ( !env )
    {
        return false;
    }

    val = env;

    if ( std::string::npos != val.find_first_of( RANK_SEPARATOR_CHAR, 0 ) )
    {
        // val = val_1,val_2,val_3

        std::vector<std::string> values;

        tokenize( values, env, RANK_SEPARATOR_CHAR );

        int pos = sRank % static_cast<int>( values.size() );

        val = values[pos];
    } 

    return true;
}

/* ----------------------------------------------------------------------------- */

void Settings::tokenize( std::vector<std::string>& values, const std::string& input, const char seperator )
{
    values.clear();

    size_t found = std::string::npos;

    do
    {
        const size_t prevFound = found + 1;
        found = input.find( seperator, prevFound );
        values.push_back( input.substr( prevFound, found - prevFound ) );
    }
    while ( found != std::string::npos );
}

bool Settings::getEnvironment( std::vector<std::string>& vals, const char* envVarName, const char separator )
{
    std::string val;

    bool found = getEnvironment( val, envVarName );

    if ( found )
    {
        tokenize( vals, val, separator );
    }

    return found;
}

/* ----------------------------------------------------------------------------- */

void Settings::printEnvironment()
{
    int i = 0;
    char *s = *environ;

    for ( ; s; i++ )
    {
        if ( strncmp( s, "SCAI_", 5 ) == 0 )
        {
            std::cout << s << std::endl;
        }

        s = *( environ + i );
    }
}

/* ----------------------------------------------------------------------------- */

void Settings::putEnvironment( const char* envVarName, const char* val, bool replace )
{
    int c_replace = 0;

    if ( replace )
    {
        c_replace = 1;
    }

    // Note: use of putenv is unsafe for auto-strings

#ifdef WIN32
    if ( replace || ( getenv( envVarName ) == NULL ) )
    {
        _putenv_s( envVarName, val );
    }
#else
    setenv ( envVarName, val, c_replace );
#endif
}

/* ----------------------------------------------------------------------------- */

void Settings::putEnvironment( const char* envVarName, int val, bool replace )
{
    std::ostringstream str_val;

    str_val << val;

    putEnvironment( envVarName, str_val.str().c_str(), replace );
}

/* ----------------------------------------------------------------------------- */

int Settings::sRank = 0;    // default value

void Settings::setRank( int rank ) 
{
    sRank = rank;
}

} /* end namespace common */

} /* end namespace scai */

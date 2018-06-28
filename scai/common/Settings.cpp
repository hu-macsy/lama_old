/**
 * @file Settings.cpp
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

extern char** environ;

namespace scai
{

namespace common
{

/* ----------------------------------------------------------------------------- */

const char* Settings::RANK_DELIMITER()
{
    return ",";
}

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

bool Settings::convertYesNoString( bool& flag, const char* stringVal )
{
    char key = static_cast<char>( toupper( stringVal[0] ) );
    bool done = true; // becomes false if no legal value has been found

    // to upper

    if ( key == '0' )
    {
        flag = false;
    }
    else if ( key == '1' )
    {
        flag = true;
    }
    else if ( key == 'J' )
    {
        flag = true;
    }
    else if ( key == 'Y' )
    {
        flag = true;
    }
    else if ( key == 'T' )
    {
        flag = true;
    }
    else if ( key == 'N' )
    {
        flag = false;
    }
    else if ( key == 'F' )
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

bool Settings::getRankedEnvironment( std::string& val, const char* envVarName )
{
    const char* env = getenv( envVarName );

    if ( !env )
    {
        return false;
    }

    val = env;

    if ( std::string::npos != val.find_first_of( RANK_DELIMITER(), 0 ) )
    {
        // val = val_1,val_2,val_3
        std::vector<std::string> values;
        tokenize( values, env, RANK_DELIMITER() );
        int pos = sRank % static_cast<int>( values.size() );
        val = values[pos];
    }

    return true;
}

/* ----------------------------------------------------------------------------- */

void Settings::tokenize( std::vector<std::string>& tokens, const std::string& input, const std::string& delimiters )
{
    tokens.clear();

    // Skip delimiters at beginning.
    std::string::size_type lastPos = input.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    std::string::size_type pos     = input.find_first_of( delimiters, lastPos );

    while ( std::string::npos != pos || std::string::npos != lastPos )
    {
        // Found a token, add it to the vector.
        tokens.push_back( input.substr( lastPos, pos - lastPos ) );
        // Skip delimiters.  Note the "not_of"
        lastPos = input.find_first_not_of( delimiters, pos );
        // Find next "non-delimiter"
        pos = input.find_first_of( delimiters, lastPos );
    }
}

bool Settings::getEnvironment( std::vector<std::string>& vals, const char* envVarName, const char* delimiters )
{
    std::string val;
    bool found = getEnvironment( val, envVarName );

    if ( found )
    {
        tokenize( vals, val, delimiters );
    }

    return found;
}

/* ----------------------------------------------------------------------------- */

void Settings::printEnvironment( std::ostream& out )
{
    int i = 0;
    char* s = *environ;

    for ( ; s; i++ )
    {
        if ( strncmp( s, "SCAI_", 5 ) == 0 )
        {
            out << s << std::endl;
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

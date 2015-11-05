/**
 * @file LibModule.cpp
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
 * @brief Implementatinos for LibModule
 *
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#include <scai/common/LibModule.hpp>

#include <scai/common/exception/Exception.hpp>

#include <dlfcn.h>
#include <dirent.h>

#include <iostream>
#include <vector>
#include <string>

#undef DEBUG_HERE

#define SUFFIX "so"

namespace scai
{

namespace common
{

static bool isDirectory( const char* dir )
{
    DIR *dp = opendir( dir );

    if ( dp == NULL )
    {
        return false;
    }

    closedir( dp );

    return true;
}

static void getFilesFromDirectory( const char* dir, const char* suffix, std::vector<std::string>& files )
{
#ifdef DEBUG_HERE
    std::cout << "getFilesfromDirectory " << dir << std::endl;
#endif

    DIR *dp = opendir( dir );

    if ( dp == NULL )
    {
        COMMON_THROWEXCEPTION( "Error (" /* << errno */ << " ) opening directory " << dir ) 
    }

    const std::string directory( dir );

    std::string slash;

    if ( *directory.rbegin() != '/' )
    {
        slash = "/";
    }
    else
    {
        slash = "";
    }

    for (;;) 
    {
        struct dirent *dirp = readdir( dp ); 
  
        if ( dirp == NULL )
        {
           break;
        }

        std::string filename = dirp->d_name;

#ifdef DEBUG_HERE
        std::cout << "File in dir " << dir << ": " << filename << std::endl;
#endif

        if ( filename.substr( filename.find_last_of( "." ) + 1 ) == suffix )
        {
            files.push_back( directory + slash + filename );
        }
    }

    closedir( dp );
}

/* -------------------------------------------------------------------------- */

LibModule::LibHandle LibModule::loadLib( const char* filename )
{
    LibHandle handle = dlopen( filename,RTLD_LAZY|RTLD_GLOBAL );

    if ( !handle )
    {
        COMMON_THROWEXCEPTION( "Cannot load library " << filename << ", " << dlerror() )
    }
    return handle;
}

/* -------------------------------------------------------------------------- */

void LibModule::freeLib( LibHandle handle )
{
    int rc = dlclose( handle );
  
    if ( rc )
    {
        COMMON_THROWEXCEPTION( "Unload library failed, " << dlerror() )
    }
};

/* -------------------------------------------------------------------------- */

void LibModule::loadLibsInDir( const char* dir )
{
    std::vector<std::string> moduleNames;

    getFilesFromDirectory( dir, SUFFIX, moduleNames );

    for ( size_t i = 0; i < moduleNames.size(); ++i )
    {

#ifdef DEBUG_HERE
        std::cout << "try to load module " << moduleNames[i] << std::endl;
#endif

        // continue if name does not match the pattern
 
        try 
        {
            loadLib( moduleNames[i].c_str() );
        }
        catch ( scai::common::Exception& e )
        {
            std::cerr << "Could not load libary " << moduleNames[i] << " in " << dir
                      << ", " << e.what() << std::endl;
        }
    }
}

/* -------------------------------------------------------------------------- */

static void tokenize( std::vector<std::string>& tokens,
                      const std::string& str,
                      const std::string& delimiters = " ")
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of( delimiters, lastPos );

    while ( std::string::npos != pos || std::string::npos != lastPos )
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr( lastPos, pos - lastPos ) );
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of( delimiters, pos );
        // Find next "non-delimiter"
        pos = str.find_first_of( delimiters, lastPos );
    }
}

/* -------------------------------------------------------------------------- */

void LibModule::loadLibsByPath( const char *path )
{
    std::vector<std::string> tokens;

    std::string input = path;
 
    tokenize( tokens, input, ":" );

    for ( size_t i = 0; i < tokens.size(); ++i )
    {
        const char* item = tokens[i].c_str();

        if ( isDirectory( item ) )
        {
            loadLibsInDir( item );
        }
        else
        {
            loadLib( item );
        }
    }
}

} /* end namespace common */

} /* end namespace scai */

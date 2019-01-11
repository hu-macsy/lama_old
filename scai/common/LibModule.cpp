/**
 * @file LibModule.cpp
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
 * @brief Implementatinos for LibModule
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#include <scai/common/LibModule.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/throw.hpp>

#include <iostream>
#include <vector>
#include <string>

#undef DEBUG_HERE

#if defined(__APPLE__)
#define SUFFIX "dylib"
#elif defined(_WIN32)
#define SUFFIX "dll"
#else // LINUX
#define SUFFIX "so"
#endif

namespace scai
{

namespace common
{

#if defined(WIN32)
// TODO: needs to be tested
static bool isDirectory( const char* dir )
{
    DWORD ftyp = GetFileAttributesA( dir );

    if ( ftyp == INVALID_FILE_ATTRIBUTES )
    {
        return false;
    }

    if ( ftyp & FILE_ATTRIBUTE_DIRECTORY )
    {
        return true;
    }

    return false;
}

/* -------------------------------------------------------------------------- */

static void getFilesFromDirectory( const char* dir, const char* suffix, std::vector<std::string>& files )
{
    WIN32_FIND_DATA ffd;
    HANDLE dp = INVALID_HANDLE_VALUE;
    dp = FindFirstFile( dir, &ffd );

    if ( dp == INVALID_HANDLE_VALUE )
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

    do
    {
        std::string filename;

        if ( ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
        {
            // todo: should we do something here?
        }
        else
        {
            filename = ffd.cFileName;
        }

#ifdef DEBUG_HERE
        std::cout << "File in dir " << dir << ": " << filename << std::endl;
#endif

        if ( filename.substr( filename.find_last_of( "." ) + 1 ) == suffix )
        {
            files.push_back( directory + slash + filename );
        }
    }
    while ( FindNextFile( dp, &ffd ) != 0 );

    FindClose( dp );
}

/* -------------------------------------------------------------------------- */

LibModule::LibHandle LibModule::loadLib( const char* filename )
{
    HINSTANCE handle = LoadLibrary( filename );

    if ( handle == NULL )
    {
        COMMON_THROWEXCEPTION( "Cannot load library " << filename ) //<< ", " << dlerror())
    }

    return handle;
}

/* -------------------------------------------------------------------------- */

void LibModule::freeLib( LibHandle handle )
{
    int rc = FreeLibrary( handle );

    if ( rc )
    {
        COMMON_THROWEXCEPTION( "Unload library failed " ) // << dlerror())
    }
}

#else
static bool isDirectory( const char* dir )
{
    DIR* dp = opendir( dir );

    if ( dp == NULL )
    {
        return false;
    }

    closedir( dp );
    return true;
}

/* -------------------------------------------------------------------------- */

static void getFilesFromDirectory( const char* dir, const char* suffix, std::vector<std::string>& files )
{
#ifdef DEBUG_HERE
    std::cout << "getFilesfromDirectory " << dir << std::endl;
#endif
    DIR* dp = opendir( dir );

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

    for ( ;; )
    {
        struct dirent* dirp = readdir( dp );

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
    LibHandle handle = dlopen( filename, RTLD_LAZY | RTLD_GLOBAL );

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

#endif

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

void LibModule::loadLibsByPath( const char* path )
{
    std::vector<std::string> tokens;
    std::string input = path;
    Settings::tokenize( tokens, input, ":" );

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

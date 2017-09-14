/**
 * @file readDir.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief readDir.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 */
/*
 * readDir.cpp
 *
 *  Created on: 31.01.2011
 *      Author: rrehrman
 */
#include <scai/benchmark/readDir.hpp>
#include <scai/benchmark/BFError.hpp>

#include <cstring>
#include <sstream>
#include <cerrno>

#if defined _WIN32 || defined _WIN64
#define __windows__ 1;
#endif

#ifndef __windows__
#include <dirent.h>
#else
#include <windows.h>
#include <strsafe.h>
#endif //__windows__

using namespace scai;

void getFilesFromDirectory( const char* dir, std::vector<std::string>& files )
{
#ifdef __unix__
    const std::string suffix = "so";
#else
#ifdef __windows__
    const std::string suffix = "dll";
#else
#ifdef __APPLE__
    const std::string suffix = "lib";
#else
#error "No platform specified. Define __unix__ for Unix, WIN32,"   \
"_WIN64 for Windows or __APPLE__ for MacOS."
#endif
#endif
#endif

#ifndef __windows__
    DIR *dp;
    struct dirent *dirp;
    if( ( dp = opendir( dir ) ) == NULL )
    {
        std::stringstream message;
        message << "Error(" << errno << ") opening " << dir;
        throw bf::BFError( message.str() );
    }

    const std::string directory( dir );
    std::string slash;

    if( directory[directory.size()] != '/' )
    {
        slash = "/";
    }
    else
    {
        slash = "";
    }

    while( ( dirp = readdir( dp ) ) != NULL )
    {
        std::string filename = dirp->d_name;
        if( filename.substr( filename.find_last_of( "." ) + 1 ) == suffix )
        {
            files.push_back( directory + slash + filename );
        }
    }
    closedir( dp );
#else
    size_t length_of_arg = 0;
    StringCchLength(dir, MAX_PATH, &length_of_arg);
    if (length_of_arg > (MAX_PATH - 3))
    {
        std::stringstream message;
        message<<"Error("<<errno<<") opening "<<dir<<" (Reason: Directory path is too long.)";
        throw bf::BFError( message.str( ) );
    }

    TCHAR szDir[MAX_PATH];

    StringCchCopy( szDir, MAX_PATH, dir );
    StringCchCat( szDir, MAX_PATH, TEXT("\\*") );

    WIN32_FIND_DATA ffd;
    HANDLE hFind = INVALID_HANDLE_VALUE;
    hFind = FindFirstFile(szDir, &ffd);

    if (INVALID_HANDLE_VALUE == hFind)
    {
        std::stringstream message;
        message<<"Error("<<errno<<") opening "<<dir<<" (Reason: FindFirstFile failed.)";
        throw bf::BFError( message.str( ) );
    }

    const std::string directory( dir );
    std::string backslash;

    if( directory[directory.size( )]!='\\' )
    {
        backslash="\\";
    }
    else
    {
        backslash="";
    }

    do
    {
        //if file is no directory
        if ( !(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) )
        {
            std::string filename = ffd.cFileName;
            if( filename.substr( filename.find_last_of( "." )+1 )==suffix )
            {
                files.push_back( directory + backslash + filename );
            }
        }
    }
    while ( FindNextFile(hFind, &ffd) != 0 );

    FindClose(hFind);

#endif //__windows__
}

void bf::getFilesFromPath( const std::string& path, std::vector<std::string>& files )
{
    char* dir;
    dir = strtok( const_cast<char*>( path.c_str() ), ":" );

    while( dir != NULL )
    {
        getFilesFromDirectory( dir, files );
        dir = strtok( NULL, ":" );
    }
}

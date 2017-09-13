/**
 * @file readDir.cpp
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
 * @brief readDir.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
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

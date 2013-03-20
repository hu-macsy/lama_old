/**
 * @file shared_lib.hpp
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
 * @brief Shared library load routines.
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */

#ifndef LAMA_SHARED_LIB_HPP_
#define LAMA_SHARED_LIB_HPP_

#ifdef WIN32
#include <Windows.h>
#define LAMA_LIB_HANDLE_TYPE HINSTANCE
#ifdef max
#undef max
#endif
#else
#include <dlfcn.h>
#define LAMA_LIB_HANDLE_TYPE void*
#endif //WIN32
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

template<typename FunctionHandleType>
int loadLibAndGetFunctionHandle( FunctionHandleType& functionHandle, LAMA_LIB_HANDLE_TYPE& handle, const char* const filename, const char* const functionName )
{
#ifdef WIN32
    handle = LoadLibrary( filename );
    if( !handle )
    {
        std::stringstream message;
        message<<"Cannot load library: "<<filename<<", because LoadLibrary failed with: "<<GetLastError();
        std::cout << message.str();
        return 1;
    }

    functionHandle = ( FunctionHandleType ) GetProcAddress ( handle , functionName );

    if( !functionHandle )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"' from '"
               <<filename<<"': "<<GetLastError();
        std::cout << message.str();
        FreeLibrary( handle );
        return 1;
    }
#else
    // load library
    handle = dlopen( filename,RTLD_LAZY|RTLD_GLOBAL );
    if( !handle )
    {
        std::stringstream message;
        message<<"Cannot load library: "<<dlerror( );
        std::cout << message.str( );
        return 1;
    }
    dlerror( );

    // getting function getInputSetRegistry( ) from loaded library.
    functionHandle = ( FunctionHandleType ) dlsym( handle,functionName );

    char* dlsym_error = NULL;
    if( ( dlsym_error=dlerror( ) )!=NULL )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"' from '"
               <<filename<<"': "<<dlsym_error;
        std::cout << message.str( );
        dlclose( handle );
        return 1;
    }
#endif //WIN32
    return 0;
}

template<typename FunctionHandleType>
int getFunctionHandle( FunctionHandleType& functionHandle, LAMA_LIB_HANDLE_TYPE& handle, const char* const functionName )
{
#ifdef WIN32
    functionHandle = ( FunctionHandleType ) GetProcAddress ( handle , functionName );

    if( !functionHandle )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"': "<<GetLastError();
        std::cout << message.str( );
        FreeLibrary( handle );
        return 1;
    }
#else
    // getting function getInputSetRegistry( ) from loaded library.
    functionHandle = ( FunctionHandleType ) dlsym( handle,functionName );

    char* dlsym_error = NULL;
    if( ( dlsym_error=dlerror( ) )!=NULL )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"': "<<dlsym_error;
        std::cout << message.str( );
        dlclose( handle );
        return 1;
    }
#endif //WIN32
    return 0;

}

inline void freeLibHandle( LAMA_LIB_HANDLE_TYPE handle )
{
#ifdef WIN32
    FreeLibrary( handle );
#else
    dlclose( handle );
#endif //WIN32
}

#endif // LAMA_SHARED_LIB_HPP_

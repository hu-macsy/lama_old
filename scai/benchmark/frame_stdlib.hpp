/**
 * @file frame_stdlib.hpp
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
 * @brief frame_stdlib.h
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

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

#include <scai/common/config.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>
#include <scai/benchmark/string_helper.hpp>

namespace scai
{

namespace bf
{

void getSharedLibraries( std::vector<std::string>& files );

template<typename FunctionHandleType,typename LibraryHandleType>
int loadLibAndGetFunctionHandle(
    FunctionHandleType& functionHandle,
    LibraryHandleType& handle,
    const char* const filename,
    const char* const functionName )
{
#ifdef WIN32
    handle = LoadLibrary( filename );
    if( !handle )
    {
        std::stringstream message;
        message<<"Cannot load library: "<<filename<<", because LoadLibrary failed with: "<<GetLastError();
        bf::BenchmarkPrinter::error( message.str( ) );
        return 1;
    }

    functionHandle = ( FunctionHandleType ) GetProcAddress ( handle , functionName );

    if( !functionHandle )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"' from '"
               <<filename<<"': "<<GetLastError();
        bf::BenchmarkPrinter::error( message.str( ) );
        FreeLibrary( handle );
        return 1;
    }
#else
    // load library
    handle = dlopen( filename, RTLD_LAZY | RTLD_GLOBAL );
    if( !handle )
    {
        std::stringstream message;
        message << "Cannot load library: " << dlerror();
        bf::BenchmarkPrinter::error( message.str() );
        return 1;
    }
    dlerror();

    // getting function getInputSetRegistry( ) from loaded library.
    functionHandle = (FunctionHandleType) dlsym( handle, functionName );

    char* dlsym_error = NULL;
    if( ( dlsym_error = dlerror() ) != NULL )
    {
        std::stringstream message;
        message << "Cannot load symbol '" << functionName << "' from '" << filename << "': " << dlsym_error;
        bf::BenchmarkPrinter::error( message.str() );
        dlclose( handle );
        return 1;
    }
#endif //WIN32
    return 0;
}

template<typename FunctionHandleType,typename LibraryHandleType>
int getFunctionHandle( FunctionHandleType& functionHandle, LibraryHandleType& handle, const char* const functionName )
{
#ifdef WIN32
    functionHandle = ( FunctionHandleType ) GetProcAddress ( handle , functionName );

    if( !functionHandle )
    {
        std::stringstream message;
        message<<"Cannot load symbol '"<<functionName<<"': "<<GetLastError();
        bf::BenchmarkPrinter::error( message.str( ) );
        FreeLibrary( handle );
        return 1;
    }
#else
    // getting function getInputSetRegistry( ) from loaded library.
    functionHandle = (FunctionHandleType) dlsym( handle, functionName );

    char* dlsym_error = NULL;
    if( ( dlsym_error = dlerror() ) != NULL )
    {
        std::stringstream message;
        message << "Cannot load symbol '" << functionName << "': " << dlsym_error;
        bf::BenchmarkPrinter::error( message.str() );
        dlclose( handle );
        return 1;
    }
#endif //WIN32
    return 0;

}

template<typename LibraryHandleType>
void freeLibHandle( LibraryHandleType handle )
{
#ifdef WIN32
    FreeLibrary( handle );
#else
    dlclose( handle );
#endif //WIN32
}

}

}

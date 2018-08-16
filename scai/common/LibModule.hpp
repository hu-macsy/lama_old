/**
 * @file LibModule.hpp
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
 * @brief Class to search and load library modules
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#pragma once

#include <scai/common/config.hpp>

#if defined( __APPLE__ )
//todo
#include <dirent.h>
#include <dlfcn.h>
#elif defined( _WIN32 )
#include <Windows.h>
#include <WinBase.h>
#include <direct.h>
#else // LINUX
#include <dlfcn.h>
#include <dirent.h>
#endif

namespace scai
{

namespace common
{

/** Static class that provides methods to find, load and unload library modules. */

class COMMON_DLL_IMPORTEXPORT LibModule
{
public:

    /** Data type definition for library handle, might be OS specific. */

#if defined( _WIN32 )
    typedef HINSTANCE LibHandle;  //!< use HINSTANCE of MSVC
#else
    typedef void* LibHandle;
#endif

    /** Load a library with its full name i.e. include suffix .so, might be absolute or relative */

    static LibHandle loadLib( const char* filename );

    /** Unload a library
     *
     *  Note: It is very likely that the libary is not unloaded now, e.g. there is no
     *        guarantee that the destructors of static objects are called.
     */
    static void freeLib( LibHandle handle );

    /** This routine reads a directory and tries to load all library modules in it that match
     *  the pattern.
     *
     *  throws an exception if directory does not exist or is not readable
     */

    static void loadLibsInDir( const char* dir );

    /** multiple directory/libaries by path string like module_1:module_2:module_directory  */

    static void loadLibsByPath( const char* path );
};

} /* end namespace common */

} /* end namespace scai */

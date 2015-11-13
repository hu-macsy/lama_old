/**
 * @file LibModule.hpp
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
 * @brief Class to search and load library modules
 *
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#pragma once

#include <dlfcn.h> 


namespace scai
{

namespace common
{

/** Static class that provides methods to find, load and unload library modules. */

class LibModule
{
public:

    /** Data type definition for library handle, might be OS specific. */

    typedef void* LibHandle;

    /** Load a library with its full name i.e. include suffix .so, might be absolute or relative */

    static LibHandle loadLib( const char* filename );

    /** Unload a library 
     *
     *  Note: It is very likely that the libary is not unloaded now, e.g. there is no
     *        guarantee that the destructors of static objects are called.
     */
    static void freeLib( LibHandle handle );

    /** This routine reads a directory and tries to load all library modules in it that match the pattern. 
     *
     *  throws an exception if directory does not exist or is not readable
     *
     */

    static void loadLibsInDir( const char* dir );

    /** multiple directory/libaries by path string like module_1:module_2:module_directory  */

    static void loadLibsByPath( const char* path );
};

} /* end namespace common */

} /* end namespace scai */

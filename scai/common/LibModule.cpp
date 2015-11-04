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

namespace scai
{

namespace common
{

LibModule::LibHandle LibModule::loadLib( const char* filename )
{
    LibHandle handle = dlopen( filename,RTLD_LAZY|RTLD_GLOBAL );

    if ( !handle )
    {
        COMMON_THROWEXCEPTION( "Cannot load library " << filename << ", " << dlerror() )
    }
    return handle;
}

void LibModule::freeLib( LibHandle handle )
{
    int rc = dlclose( handle );
  
    if ( rc )
    {
        COMMON_THROWEXCEPTION( "Unload library failed, " << dlerror() )
    }
};

} /* end namespace common */

} /* end namespace scai */

/**
 * @file examples/UseModule.cpp
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
 * @brief ToDo: Missing description in ./examples/UseModule.cpp
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#include "DynRoutine.hpp"

#include <scai/common/LibModule.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <memory>

using std::unique_ptr;
using namespace scai::common;

void runIt()
{
    // as long as module is not loaded,  Function1 is not registered at factory
    if ( !DynRoutine::canCreate( "Function1" ) )
    {
        std::cout << "Function1 not available" << std::endl;
        return;
    }

    std::cout << "Function1 can be created." << std::endl;
    unique_ptr<DynRoutine> d( DynRoutine::create( "Function1" ) );
    d->doIt();
}

int main( int argc, const char** argv )
{
    runIt();  // shoud not work
    std::string loadPath;

    if ( Settings::getEnvironment( loadPath, "SCAI_LIBRARY_PATH" ) )
    {
        std::cout << "Load all module libraries in loadPath " << loadPath << std::endl;
        LibModule::loadLibsByPath( loadPath.c_str() );
        runIt();   // should now work
    }
    else
    {
        if ( argc != 2 )
        {
            COMMON_THROWEXCEPTION( "UseModule <libname>" )
        }

        const char* lib = argv[1];

        std::cout << "Load " << lib << std::endl;

        LibModule::LibHandle handle = LibModule::loadLib( lib );

        runIt();

        std::cout << "Library loaded successfully, now unload" << std::endl;

        LibModule::freeLib( handle );

        std::cout << "Library freed, still try to use it" << std::endl;

        runIt();   // undefined whether it works
    }
}

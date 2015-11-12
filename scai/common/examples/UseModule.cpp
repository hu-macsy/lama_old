
#include "DynRoutine.hpp"

#include <scai/common/LibModule.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

#include <iostream>

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

    unique_ptr<DynRoutine> d = DynRoutine::create( "Function1" );

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

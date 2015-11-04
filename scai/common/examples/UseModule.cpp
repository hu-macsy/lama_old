
#include <scai/common/LibModule.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/common/Settings.hpp>

#include "DynRoutine.hpp"

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

    DynRoutine* d = DynRoutine::create( "Function1" );

    d->doIt();

    delete d;
}

int main( int, const char** )
{
    runIt();  // shoud not work

    std::string directory;

    if ( Settings::getEnvironment( directory, "SCAI_LIBRARY_PATH" ) )
    {
        std::cout << "Load all module libraries in directory " << directory << std::endl;

        LibModule::loadLibsInDir( directory.c_str(), "" );

        runIt();   // should now work
    }
    else
    {
        std::cout << "Load libmodule.s0" << std::endl;

        LibModule::LibHandle handle = LibModule::loadLib( "libmodule.so" );

        runIt();   // should now work

        LibModule::freeLib( handle );

        std::cout << "Library freed, still try to use it" << std::endl;

        runIt();   // undefined whether it works
    }
}

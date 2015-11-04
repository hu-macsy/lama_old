
#include <scai/common/LibModule.hpp>
#include <scai/common/exception/Exception.hpp>

#include "DynRoutine.hpp"

using namespace scai::common;

int main( int, const char** )
{
    // this should not be available yet

    std::cout << "Function1 available (before loadLib): " << DynRoutine::canCreate( "Function1" ) << std::endl;

    LibModule::LibHandle handle = LibModule::loadLib( "libmodule.so" );

    std::cout << "Function1 available (after loadLib): " << DynRoutine::canCreate( "Function1" ) << std::endl;

    // Now we can get it

    DynRoutine* d = DynRoutine::create( "Function1" );

    // We can use the object created by routines of the library module

    d->doIt();

    LibModule::freeLib( handle );

    // unclear what happens because entry in factory has not been deleted

    std::cout << "Function1 available (freed loadLib): " << DynRoutine::canCreate( "Function1" ) << std::endl;
}


#include "DynRoutine.hpp"

#include <iostream>

using namespace scai::common;

/** Define a new class that will be added to the factory for DynRoutine. */

class Function1 : public DynRoutine, public DynRoutine::Register<Function1>
{
public:

    /** This is the implementation for the virtual function of base class DynRoutine */

    void doIt()
    {
        std::cout << "Function1: doit" << std::endl;
    }

    /** This static functions provides the value for which this class registers in the factory. */

    static inline std::string createValue()
    {   
        return "Function1";
    }
    
    /** Method that creates objects of type Function1 that will be used for registration. */

    static DynRoutine* create()
    {
        return new Function1();
    }
};

DynRoutine::Register<Function1>::RegisterGuard DynRoutine::Register<Function1>::registerGuard;

/** Guard class to call functions for load/unload of the module */

struct CGuard
{
    CGuard()
    {
        std::cout << "Module is loaded" << std::endl;
    }

    ~CGuard()
    {
        std::cout << "Module is unloaded" << std::endl;
    }
};

static CGuard ModuleGuard;

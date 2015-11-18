/**
 * @file common/examples/DummyModule.cpp
 *
 * @brief Example of a very simple library module
 *
 * @author Thomas Brandes
 * @date 19.10.2015
 */

#include <iostream>

// own namespace avoids conflict with any other libraries

namespace dummy_module

{

/** Guard class to call functions for load/unload of the module */

struct CGuard
{
    CGuard()
    {
        std::cout << "DummyModule is loaded" << std::endl;
    }

    ~CGuard()
    {
        std::cout << "DummyModule is unloaded" << std::endl;
    }
};

static CGuard Module1Guard;

}


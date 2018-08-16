/**
 * @file examples/Module.cpp
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
 * @brief ToDo: Missing description in ./examples/Module.cpp
 * @author Thomas Brandes
 * @date 04.11.2015
 */

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

template DynRoutine::Register<Function1>::RegisterGuard DynRoutine::Register<Function1>::registerGuard;

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

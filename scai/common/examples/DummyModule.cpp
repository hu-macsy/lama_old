/**
 * @file common/examples/DummyModule.cpp
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
 * @brief Example of a very simple library module
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


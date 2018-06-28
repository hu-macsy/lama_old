/**
 * @file common/examples/DemoMetaProgramming.cpp
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
 * @brief Examples for using Template Type lists
 * @author Thomas Brandes
 * @date 15.03.2016
 */

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

#include <iostream>
#include <typeinfo>

// define the list of types

#define MY_TYPELIST SCAI_TYPELIST( SCAI_NUMERIC_TYPES_HOST )

using namespace scai;
using namespace common;

// This method is instantiated and called for all types of MY_TYPELIST

template<typename T>
void output()
{
    bool isCmp = isComplex( TypeTraits<T>::stype );
    std::cout << "Running output<" << TypeTraits<T>::id() << ">, is complex = " << isCmp << std::endl;
    T x = 1;
    x = x / 3;
    std::cout << TypeTraits<T>::stype << " : " << x << std::endl;
}

// general defintion of Calling for template list

template<typename TList> struct Calling;

// termination call

template<> struct Calling<mepr::NullType>
{
    static void call()
    {
    }
};

// call output for header T and recursive call for tail of list

template<typename HeadType, typename TailTypes>
struct Calling<mepr::TypeList<HeadType, TailTypes> >
{
    static void call()
    {
        output<HeadType>();
        Calling<TailTypes>::call();
    }
};

int main()
{
    Calling<MY_TYPELIST>::call();
}

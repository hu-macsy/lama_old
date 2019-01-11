/**
 * @file common/examples/DemoFactory1.cpp
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
 * @brief Example of a factory
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Factory1.hpp>
#include <scai/common/Printable.hpp>

#include <iostream>

using namespace scai::common;

enum Kind
{
    D1,
    D2,
    D3
};

class Base : public Factory1<Kind, int, Base*>, public Printable
{
public:
    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "Base";
    }
};

class Derived1 : public Base //  Factory<int, *Base>::Manager<Derived1, 1>
{
public:

    int mVal;

    Derived1( int val ) : mVal( val ) {}

    static Base* create( int val )
    {
        return new Derived1( val );
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "Derived1, mVal = " << mVal;
    }

    static bool init();
    static bool initialized;

};

bool Derived1::initialized = Derived1::init();

bool Derived1::init()
{
    addCreator( D1, &create );
    return true;
}

class Derived2 : public Base, private Factory1<Kind, int, Base*>::Register<Derived2>
{
public:

    int mVal;

    Derived2( int val ) : mVal( val ) {}

    static inline Kind createValue()
    {
        return D2;
    }

    static Base* create( int val )
    {
        return new Derived2( val );
    }

private:

    virtual void writeAt( std::ostream&  stream ) const
    {
        stream << "Derived2, mVal = " << mVal;
    }
};

// Some compilers require explicit instantiation of the register guard

template Base::Register<Derived2>::RegisterGuard Base::Register<Derived2>::registerGuard;

int main()
{
    std::vector<Kind> values;
    Base::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::cout << "Registered values[" << i << "] = " << values[i] << std::endl;
    }

    Base* obj1 = Base::create( D1, 15 );
    Base* obj2 = Base::create( D2, -5 );
    std::cout << "obj1 is " << *obj1 << std::endl;
    std::cout << "obj2 is " << *obj2 << std::endl;
}


/**
 * @file common/examples/DemoFactory1.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Example of a factory
 *
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <common/Factory.hpp>
#include <common/Printable.hpp>

#include <iostream>

using namespace common;

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
        return new Derived2( val);
    }

private:

    virtual void writeAt( std::ostream&  stream ) const
    {
        stream << "Derived2, mVal = " << mVal;
    }
};

int main()
{
    std::vector<Kind> values;

    Base::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::cout << "Registered values[" << i << "] = " << values[i] << std::endl;
    }

    Base* obj1 = Base::create(D1, 15);
    Base* obj2 = Base::create(D2, -5);

    std::cout << "obj1 is " << *obj1 << std::endl;
    std::cout << "obj2 is " << *obj2 << std::endl;
}


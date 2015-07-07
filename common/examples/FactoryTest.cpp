/**
 * @file common/examples/FactoryTest.cpp
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
    D2
};

class Base : public Factory<Kind, Base*>, public Printable
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

    static Base* create()
    {
        return new Derived1();
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "Derived1";
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

class Derived2 : public Base, private Factory<Kind, Base*>::Register<Derived2, D2>
{
public:
    static Base* create()
    {
        return new Derived2();
    }

private:

    virtual void writeAt( std::ostream&  stream ) const
    {
        stream << "Derived2";
    }
};

int main()
{
    Base* obj1 = Base::create(D1);
    Base* obj2 = Base::create(D2);

    std::cout << "obj1 is " << *obj1 << std::endl;
    std::cout << "obj2 is " << *obj2 << std::endl;
}


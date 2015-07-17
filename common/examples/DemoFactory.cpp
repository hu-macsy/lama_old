/**
 * @file common/examples/DemoFactory.cpp
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
#include <typeinfo>

using namespace common;

/** Base class that provides by deriving from Factory a factory with a create routine.
 *
 *  The input value for the create routine is a string, the output value a pointer to
 *  a Base object.
 *
 *  \code
 *       *Base create( std::string )
 *  \endcode
 */

class Base : public Factory<std::string, Base*>, public Printable
{
public:

    // Override default method of class Printable

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "Base";
    }
};

/** Derived class that registers manually in the factory. */

class Derived1 : public Base //  Factory<int, *Base>::Manager<Derived1, 1>
{
public:

    /** Method that creates objects of type Derived2 that will be used for registration. */

    static Base* create()
    {
        return new Derived1();
    }

private:

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
    addCreator( "D1", &create );
    return true;
}

/** @brief Derived class of Base that uses Register class of Factory.
 *
 *  This class demonstrates automatic registration in the factory of the base class.
 */
class Derived2 : public Base, Base::Register<Derived2>
{
public:

    /** This static functions provides the value for which this class registers in the factory. */

    static inline std::string createValue()
    {
        return "D2";
    }

    /** Method that creates objects of type Derived2 that will be used for registration. */

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

/** @brief Derived template class of Base that uses Register class of Factory.
 */
template<typename T>
class Derived : public Base, Base::Register<Derived<T> >
{
public:

    /** This static functions provides the value for which this class registers in the factory. */

    static inline std::string createValue()
    {
        return typeid( T ).name();
    }

    /** Method that creates objects of type Derived2 that will be used for registration. */

    static Base* create()
    {
        return new Derived<T>();
    }

private:

    T mValue;

    virtual void writeAt( std::ostream&  stream ) const
    {
        stream << "Derived<" << typeid( T ).name() << ">";
    }
};

template class Derived<int>;
template class Derived<float>;

using namespace std;

int main()
{
    vector<string> values;  // string is create type for the factory

    Base::getCreateValues( values );

    cout << "Factory of Base: " << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "   Registered values[" << i << "] = " << values[i] << endl;
    }

    Base* obj1 = Base::create( "D1" );
    Base* obj2 = Base::create( "D2" );
    Base* obj3 = Base::create( typeid( int ).name() );

    cout << "obj1 is " << *obj1 << endl;
    cout << "obj2 is " << *obj2 << endl;

    cout << "obj3 is " << *obj3 << endl;

    try 
    {
        Base* obj = Base::create( "na" );
        cout << "obj is " << *obj << endl;
    }
    catch ( Exception& ex )
    {
        cout << "Caught exception: " << ex.what() << endl;
    }
}


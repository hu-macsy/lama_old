/**
 * @file common/examples/DemoFactory.cpp
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

#include <scai/common/Factory.hpp>
#include <scai/common/Printable.hpp>

#include <iostream>
#include <typeinfo>

using namespace scai::common;

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

/* -----------------------------------------------------------------------------*/

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

/* -----------------------------------------------------------------------------*/


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

// Explicit instantiation of static CRTP variabled required for some compilers
// Does not work on all compilers: template class Base::Register<Derived2>;

template Base::Register<Derived2>::RegisterGuard Base::Register<Derived2>::registerGuard;

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

// Register guard instantiation, does also template instantiation

template Base::Register<Derived<int> >::RegisterGuard
Base::Register<Derived<int> >::registerGuard;
template Base::Register<Derived<float> >::RegisterGuard
Base::Register<Derived<float> >::registerGuard;

using namespace std;

/* -----------------------------------------------------------------------------*/

int main()
{
    vector<string> values;  // string is create type for the factory
    Base::getCreateValues( values );
    cout << "Factory of Base: " << values.size() << " entries" << endl;

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

    cout << "DemoFactory terminated correctly" << endl;
}


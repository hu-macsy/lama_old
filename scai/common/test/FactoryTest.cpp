/**
 * @file FactoryTest.cpp
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
 * @brief Test routines for class Factory
 * @author Thomas Brandes
 * @date 05.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Factory.hpp>
#include <scai/common/Printable.hpp>

using namespace scai::common;

/** Base class that provides by deriving from Factory a factory with a create routine.
 *  As it is also derived from Printable, it becomes polymorphic and can be used
 *  with the dynamic cast operator.
 */

class Base : public Factory<std::string, Base*>, public Printable
{
};

/* -----------------------------------------------------------------------------*/

/** @brief Derived class of Base that uses Register class of Factory.  */

class Derived : public Base, Base::Register<Derived>
{
public:

    static inline std::string createValue()
    {
        return "D";
    }

    /** Method that creates objects of type Derived that will be used for registration. */

    static Base* create()
    {
        return new Derived();
    }
};

// register variable

template Base::Register<Derived>::RegisterGuard Base::Register<Derived>::registerGuard;

/* -----------------------------------------------------------------------------*/

/** @brief Derived template class of Base that uses Register class of Factory.  */

template<typename T>
class TDerived : public Base, Base::Register<TDerived<T> >
{
public:

    /** This static functions provides the value for which this class registers in the factory. */

    static inline std::string createValue()
    {
        return typeid( T ).name();
    }

    /** Method that creates objects of type TDerived that will be used for registration. */

    static Base* create()
    {
        return new TDerived<T>();
    }
};

// Register guard instantiation, does also template instantiation

template Base::Register<TDerived<int> >::RegisterGuard Base::Register<TDerived<int> >::registerGuard;
template Base::Register<TDerived<float> >::RegisterGuard Base::Register<TDerived<float> >::registerGuard;

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE( FactoryTest )
{
    std::vector<std::string> values;  // string is create type for the factory
    Base::getCreateValues( values );
    size_t size_expected = 3;  // we have registered 3 classes in factory
    BOOST_CHECK_EQUAL( size_expected, values.size() );
    BOOST_CHECK( Base::canCreate( "D" ) );
    BOOST_CHECK( !Base::canCreate( "F" ) );
    BOOST_CHECK( Base::canCreate( typeid( int ).name() ) );
    BOOST_CHECK( !Base::canCreate( typeid( double ).name() ) );
    BOOST_CHECK_THROW( { Base::create( "e" ); }, Exception );
    Base* obj = Base::create( "D" );
    Derived* derivedObj = dynamic_cast<Derived*>( obj );
    BOOST_CHECK( derivedObj );
    obj = Base::create( typeid( int ).name() );
    TDerived<int>* intObj = dynamic_cast<TDerived<int>*>( obj );
    TDerived<float>* floatObj = dynamic_cast<TDerived<float>*>( obj );
    TDerived<double>* doubleObj = dynamic_cast<TDerived<double>*>( obj );
    BOOST_CHECK( intObj );
    BOOST_CHECK( floatObj == NULL );
    BOOST_CHECK( doubleObj == NULL );
}

/**
 * @file Factory1Test.cpp
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
 * @brief Test routines for class Factory1
 *
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Factory1.hpp>
#include <scai/common/Printable.hpp>

using namespace scai::common;

/** Base class that provides by deriving from Factory a factory with a create routine. 
 *  As it is also derived from Printable, it becomes polymorphic and can be used
 *  with the dynamic cast operator.
 */

class Base : public Factory1<std::string, int, Base*>, public Printable
{
};

/* -----------------------------------------------------------------------------*/

/** @brief Derived class of Base that uses Register class of Factory1.  */

class Derived : public Base, Base::Register<Derived>
{
public:

    /** Constructor uses InputType == int */

    Derived( int val ) : mValue( val )
    {
    }

    static inline std::string createValue()
    {
        return "D";
    }

    /** Method that creates objects of type Derived that will be used for registration. */

    static Base* create( int val )
    {
        return new Derived( val );
    }

    int getVal() const 
    {
        return mValue;
    }

private:

    int mValue;
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
    
    static Base* create( int )
    {   
        return new TDerived<T>();
    }
};

// Register guard instantiation, does also template instantiation
    
template Base::Register<TDerived<int> >::RegisterGuard Base::Register<TDerived<int> >::registerGuard;
template Base::Register<TDerived<float> >::RegisterGuard Base::Register<TDerived<float> >::registerGuard;

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE( Factory1Test )
{
    std::vector<std::string> values;  // string is create type for the factory

    Base::getCreateValues( values );

    size_t size_expected = 3;  // we have registered 3 classes in factory

    BOOST_CHECK_EQUAL( size_expected, values.size() );

    BOOST_CHECK( Base::canCreate( "D" ) );
    BOOST_CHECK( !Base::canCreate( "F" ) );
    BOOST_CHECK( Base::canCreate( typeid(int).name() ) );
    BOOST_CHECK( !Base::canCreate( typeid(double).name() ) );

    BOOST_CHECK_THROW( { Base::create( "e", 1 ); }, Exception );

    Base* obj = Base::create( "D", 5 );
    Derived* derivedObj = dynamic_cast<Derived*>( obj );
    BOOST_REQUIRE( derivedObj );
    BOOST_CHECK_EQUAL( 5, derivedObj->getVal() );

    obj = Base::create( typeid(int).name(), 0 );
    TDerived<int>* intObj = dynamic_cast<TDerived<int>*>( obj );
    TDerived<float>* floatObj = dynamic_cast<TDerived<float>*>( obj );
    TDerived<double>* doubleObj = dynamic_cast<TDerived<double>*>( obj );

    BOOST_CHECK( intObj );
    BOOST_CHECK( floatObj == NULL );
    BOOST_CHECK( doubleObj == NULL );
}

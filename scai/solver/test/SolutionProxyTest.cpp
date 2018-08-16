/**
 * @file SolutionProxyTest.cpp
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
 * @brief Test routines for the class SolutionProxy.
 * @author Matthias Makulla, Jiri Kraus
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/solver/SolutionProxy.hpp>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;

typedef DefaultReal ValueType;
typedef Vector<ValueType> VectorType;
typedef SolutionProxy<ValueType> ProxyType;

/* --------------------------------------------------------------------- */

/** Test configuration generated for each individual test. */

struct SolutionProxyTestConfig
{
    SolutionProxyTestConfig() : 

        vPtr( new DenseVector<ValueType>() ),
        mProxy( vPtr.get() )
    {
        vPtr->setSameValue( 3, -5 );
    }

    ~SolutionProxyTestConfig()
    {
    }

    std::unique_ptr<DenseVector<ValueType> > vPtr;
    ProxyType mProxy;
};

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( SolutionProxyTest , SolutionProxyTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.SolutionProxyTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testOperators )
{
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    mProxy.setDirty( false );
    BOOST_CHECK( ! mProxy.isDirty() );
    VectorType& vec = mProxy.getReference();
    vec[0] = 3;
    // using reference makes proxy dirty
    BOOST_CHECK( mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testSetAndIsDirty )
{
    mProxy.setDirty( true );
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    mProxy.setDirty( false );
    BOOST_CHECK_EQUAL( false, mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testGetConstReference )
{
    mProxy.setDirty( false );
    const VectorType& vec = mProxy.getConstReference();
    BOOST_CHECK( !mProxy.isDirty() );
    ValueType s = vec[0];
    BOOST_CHECK_EQUAL( ValueType( -5 ), s );
    // using const reference makes proxy not dirty
    BOOST_CHECK( ! mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

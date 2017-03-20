/**
 * @file SolutionProxyTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Test routines for the class SolutionProxy.
 * @author Matthias Makulla, Jiri Kraus
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/solver/SolutionProxy.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

typedef Vector VectorType;
typedef SolutionProxy ProxyType;

/* --------------------------------------------------------------------- */

struct SolutionProxyTestConfig
{
    SolutionProxyTestConfig() : mProxy( new DenseVector<RealType>( 3, -5.0 ) )
    {
    }

    ~SolutionProxyTestConfig()
    {
        // the vector must be deleted explicitly
        VectorType& vec = mProxy.getReference();
        delete &vec;
    }

    ProxyType mProxy;
};

BOOST_FIXTURE_TEST_SUITE( SolutionProxyTest , SolutionProxyTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.SolutionProxyTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testOperators )
{
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    mProxy.setDirty( false );
    BOOST_CHECK( ! mProxy.isDirty() );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), (*mProxy).getValue( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), (*mProxy).getValue( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), (*mProxy).getValue( 2 ) );
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
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 2 ) );
    // using const reference makes proxy not dirty
    BOOST_CHECK( ! mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testGetReference )
{
    mProxy.setDirty( false );
    BOOST_CHECK( !mProxy.isDirty() );
    VectorType& vec = mProxy.getReference();
    BOOST_CHECK( mProxy.isDirty() );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec.getValue( 2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

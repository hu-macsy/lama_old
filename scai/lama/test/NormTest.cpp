/**
 * @file NormTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Contains general tests for each Norm class registered in Norm factory.
 * @author Thomas Brandes
 * @date 16.06.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/Norm.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/common/test/TestMacros.hpp>

typedef SCAI_TEST_TYPE ValueType;

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

/** Class for vector of all registered Norm objects in factory */

class Norms : public std::vector<lama::NormPtr>
{

public:

    /** Constructor creates already the list with all storage pointers. */

    Norms()
    {
        using namespace scai::lama;

        std::vector<std::string> values;  //  all create values

        Norm::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            NormPtr normPtr( scai::lama::Norm::create( values[i] ) );
            push_back( normPtr );
        }
    }

    // Destructor will free all matrix storages due to use of shared pointers
};

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( NormTest );

/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_LOGGER( logger, "Test.NormTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( positiveHomogeneityTest )
{
    scai::lama::DenseVector<ValueType> x( 4, 1.0 );
    scai::lama::Scalar s = 3.0;

    scai::lama::DenseVector<ValueType> tmp( s* x );

    Norms allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm& norm = *allNorms[i];

        // Homogeneity test
 
        BOOST_CHECK_EQUAL( norm.apply ( tmp ), s * norm.apply( x ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( triangleInequalityTest )
{
    scai::lama::DenseVector<ValueType> x( 2, 2.0 );
    scai::lama::DenseVector<ValueType> y( 2, 2.0 );
    scai::lama::DenseVector<ValueType> z( x + y );

    Norms allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm& norm = *allNorms[i];

        // Inequality test
 
        Scalar nz = norm.apply( z );
        Scalar nxy = norm.apply( x ) + norm.apply( y );

        BOOST_CHECK( nz.getValue<ValueType>() <= nxy.getValue<ValueType>() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( zeroVectorTest )
{
    scai::lama::DenseVector<ValueType> x( 4, 0.0 );

    Norms allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm& norm = *allNorms[i];

        // zero test

        BOOST_CHECK_EQUAL( norm.apply( x ), Scalar( 0 ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( zeroMatrixTest )
{
    scai::lama::DenseMatrix<ValueType> m( 3, 4 );

    m.getLocalStorage().setZero();

    Norms allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm& norm = *allNorms[i];

        // zero test

        BOOST_CHECK_EQUAL( norm.apply( m ), Scalar( 0 ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( writeTest )
{
    Norms allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm& norm = *allNorms[i];

        std::ostringstream out;
        std::ostringstream outB;
        std::ostringstream outD;

        out << norm;
        norm.Norm::writeAt( outB );
        norm.writeAt( outD );

        BOOST_CHECK( out.str().length() > 0 );
        BOOST_CHECK( outD.str() == out.str() );
        BOOST_CHECK( outB.str() != outD.str() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();


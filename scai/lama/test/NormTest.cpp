/**
 * @file NormTest.cpp
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
 * @brief Contains general tests for each Norm class registered in Norm factory.
 * @author Thomas Brandes
 * @date 16.06.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/norm/Norm.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/common/test/TestMacros.hpp>

typedef SCAI_TEST_TYPE ValueType;

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

/** Class for vector of all registered Norm objects in factory */

template<typename ValueType>
class Norms : public std::vector<lama::NormPtr<ValueType> >
{

public:

    /** Constructor creates already the list with all storage pointers. */

    Norms()
    {
        using namespace scai::lama;

        std::vector<std::string> values;  //  all create values

        Norm<ValueType>::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            NormPtr<ValueType> normPtr( Norm<ValueType>::create( values[i] ) );
            this->push_back( normPtr );
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
    auto x = fill<DenseVector<ValueType>>( 4, 1.0 );

    ValueType s = 3;
    auto tmp = eval<DenseVector<ValueType>>( s * x );

    Norms<ValueType> allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm<ValueType>& norm = *allNorms[i];

        // Homogeneity test

        BOOST_CHECK_EQUAL( norm.apply ( tmp ), s * norm.apply( x ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( triangleInequalityTest )
{
    auto x = fill<DenseVector<ValueType>>( 2, 2.0 );
    auto y = fill<DenseVector<ValueType>>( 2, 2.0 );
    auto z = eval<DenseVector<ValueType>>( x + y );

    Norms<ValueType> allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm<ValueType>& norm = *allNorms[i];

        // Inequality test

        auto nz  = norm.apply( z );
        auto nxy = norm.apply( x ) + norm.apply( y );

        BOOST_CHECK( nz <= nxy );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( zeroVectorTest )
{
    auto x = fill<DenseVector<ValueType>>( 4, 0.0 );

    Norms<ValueType> allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm<ValueType>& norm = *allNorms[i];

        // zero test

        BOOST_CHECK_EQUAL( norm.apply( x ), 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( zeroMatrixTest )
{
    auto m = zero<DenseMatrix<ValueType>>( 3, 4 );

    Norms<ValueType> allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm<ValueType>& norm = *allNorms[i];

        // zero test

        BOOST_CHECK_EQUAL( norm.apply( m ), 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( writeTest )
{
    Norms<ValueType> allNorms;

    for ( size_t i = 0; i < allNorms.size(); ++i )
    {
        Norm<ValueType>& norm = *allNorms[i];

        std::ostringstream out;
        std::ostringstream outB;
        std::ostringstream outD;

        out << norm;
        norm.Norm<ValueType>::writeAt( outB );
        norm.writeAt( outD );

        BOOST_CHECK( out.str().length() > 0 );
        BOOST_CHECK( outD.str() == out.str() );
        BOOST_CHECK( outB.str() != outD.str() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();


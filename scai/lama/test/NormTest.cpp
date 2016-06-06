/**
 * @file NormTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Contains the implementation of the class NormTest
 * @author Alexander Büchel, Micha
 * @date 03.02.2012
 */

#include <scai/lama/test/NormTest.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

typedef SCAI_TEST_TYPE ValueType;

LAMA_COMMON_TEST_CASE( NormTest, positiveHomogeneityTest )
scai::lama::DenseVector<ValueType> x( 4, 1.0 );
scai::lama::Scalar s = 3.0;

scai::lama::DenseVector<ValueType> tmp( s* x );

//Inequality test
BOOST_CHECK_EQUAL( mNorm.apply ( tmp ), s* mNorm.apply( x ) );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE( NormTest, triangleInequalityTest )
scai::lama::DenseVector<ValueType> x( 2, 2.0 );
scai::lama::DenseVector<ValueType> y( 2, 2.0 );
scai::lama::DenseVector<ValueType> z( x + y );

BOOST_CHECK( mNorm.apply( z ) == mNorm.apply( x ) + mNorm.apply( y ) );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE( NormTest, ZeroVectorTest )
scai::lama::DenseVector<ValueType> x( 4, 0.0 );
BOOST_CHECK_EQUAL( mNorm.apply( x ), 0.0 );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE_RUNNER( NormTest )
{
    positiveHomogeneityTest();
    triangleInequalityTest();
    ZeroVectorTest();
}

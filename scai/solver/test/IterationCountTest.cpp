/**
 * @file IterationCountTest.cpp
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
 * @brief Contains the implementation of the class IterationCountTest.
 * @author Matthias Makulla, Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( IterationCountTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.IterationCountTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    IterationCount<ValueType> criterion1( 5 );
    IterationCount<ValueType> criterion2( 10 );
    IterationCount<ValueType> criterion3;
    IterationCount<ValueType> criterion4( criterion2 );

    BOOST_CHECK_EQUAL( criterion1.getIterationExtrema(), IndexType( 5 ) );
    BOOST_CHECK_EQUAL( criterion2.getIterationExtrema(), IndexType( 10 ) );
    BOOST_CHECK_EQUAL( criterion3.getIterationExtrema(), IndexType( 1 ) );
    BOOST_CHECK_EQUAL( criterion4.getIterationExtrema(), criterion2.getIterationExtrema() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetAndGetIterationCountTest, ValueType, scai_numeric_test_types )
{
    std::unique_ptr<IterationCount<ValueType> > criterion1( new IterationCount<ValueType>( 5 ) );
    std::unique_ptr<IterationCount<ValueType> > criterion2( new IterationCount<ValueType>( 10 ) );

    criterion1->setIterationExtrema( 15 );
    criterion2->setIterationExtrema( 3 );

    BOOST_CHECK_EQUAL( IndexType( 15 ), criterion1->getIterationExtrema() );
    BOOST_CHECK_EQUAL( IndexType( 3 ), criterion2->getIterationExtrema() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testIsSatisfied, ValueType, scai_numeric_test_types )
{
    const IndexType maxIterations = 10;

    Jacobi<ValueType> jacobi( "IterationCountTest Jacobi" );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();

    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );

    auto solution = denseVectorFill<ValueType>( 3, 1.0 );

    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( maxIterations ) );
    jacobi.initialize( coefficients );
    jacobi.setStoppingCriterion( criterion );
    jacobi.solve( solution, rhs );
    BOOST_CHECK_EQUAL( maxIterations, jacobi.getIterationCount() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( writeAtTest, ValueType, scai_numeric_test_types )
{
    IterationCount<ValueType> criterion( 5 );
    SCAI_COMMON_WRITEAT_TEST( criterion );
    std::stringstream mStream;
    mStream << criterion;
    std::string mString = mStream.str();
    BOOST_CHECK( mString.compare( "Maximal" ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

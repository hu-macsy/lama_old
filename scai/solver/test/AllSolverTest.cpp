/**
 * @file AllSolverTest.cpp
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
 * @brief Generic tests applied for each solver class provided by the solver factory.
 * @author Thomas Brandes
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/solver/Solver.hpp>
#include <scai/solver/test/TestMacros.hpp>

using scai::IndexType;

using namespace scai::lama;
using namespace scai::solver;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( AllSolverTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.AllSolverTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( writeAtTest, ValueType, scai_numeric_test_types )
{
    // Get all available solvers

    std::vector<std::string> values;

    Solver<ValueType>::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); i++ )
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr<ValueType> solver( Solver<ValueType>::getSolver( values[i] ) );
        SCAI_COMMON_WRITEAT_TEST( *solver );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( copyTest, ValueType, scai_numeric_test_types )
{
    // Get all available solvers

    std::vector<std::string> values;
    Solver<ValueType>::getCreateValues( values );
    const int numSolvers = ( int )values.size();

    for ( int i = 0; i < numSolvers; i++ )
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr<ValueType> solver1( Solver<ValueType>::getSolver( values[i] ) );
        SolverPtr<ValueType> solver2( solver1->copy() );
        BOOST_CHECK_EQUAL( solver1->getId(), solver2->getId() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( solveWithoutInitialization, ValueType, scai_numeric_test_types )
{
    // Some test data are required to call the solve method, we don't want to work with these data
    // therefore we do not need to set a context or a proper distribution here

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator::buildPoisson2D( coefficients, 9, N1, N2 );

    const ValueType solutionInitValue = 1.0;
    auto colDist = coefficients.getColDistributionPtr();

    auto solution      = denseVectorFill<ValueType>( colDist, solutionInitValue );
    auto exactSolution = denseVectorFill<ValueType>( colDist, solutionInitValue + 1 );
    auto rhs           = denseVectorEval( coefficients * exactSolution );

    // Test for all registered solvers ( with this value type )

    std::vector<std::string> values;
    Solver<ValueType>::getCreateValues( values );
    const int numSolvers = ( int )values.size();

    for ( int i = 0; i < numSolvers; i++ )
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr<ValueType> solver( Solver<ValueType>::getSolver( values[i] ) );
        BOOST_CHECK_THROW ( {solver->solve( solution, rhs );}, scai::common::Exception );
    }
}

BOOST_AUTO_TEST_SUITE_END();

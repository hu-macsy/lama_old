/**
 * @file SolverTest1.cpp
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
 * @brief SolverTest1.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/solver/Solver.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::solver;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SolverTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SolverTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr solver( Solver::create( values[i], "" ) );

        SCAI_COMMON_WRITEAT_TEST( *solver );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr solver( Solver::create( values[i], "" ) );

        SolverPtr solverCpyPtr = solver->copy();
        // TODO: do a more proper test here!
        BOOST_CHECK_EQUAL( solverCpyPtr->getId(), "" );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( solveWithoutInitialization )
{
    typedef SCAI_TEST_TYPE ValueType;

    // Some test data are required to call the solve method, we don't want to work with these data
    // therefore we do not need to set a context or a proper distribution here
    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );

    const ValueType solutionInitValue = 1.0;
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), solutionInitValue );

    DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), solutionInitValue+1.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        SolverPtr solver( Solver::create( values[i], "" ) );

        BOOST_CHECK_THROW ( {solver->solve( solution, rhs );}, scai::common::Exception );
    }
}

BOOST_AUTO_TEST_SUITE_END();

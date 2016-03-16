/**
 * @file SolverTest1.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief SolverTest1.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 * @since 2.0.0
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
    for(int i=0; i < numSolvers; i++){
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        Solver* solver = Solver::create( values[i], "" );

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
    for(int i=0; i < numSolvers; i++){
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        Solver* solver = Solver::create( values[i], "" );

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
    for(int i=0; i < numSolvers; i++){
        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );
        Solver* solver = Solver::create( values[i], "" );

        BOOST_CHECK_THROW ( {solver->solve( solution, rhs );}, scai::common::Exception );
    }
}

BOOST_AUTO_TEST_SUITE_END();

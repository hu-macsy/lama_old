/**
 * @file IterativeSolverTest.cpp
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
 * @endlicense
 *
 * @brief IterativeSolverTest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/IterativeSolver.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace scai::lama;
using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( IterativeSolverTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.IterativeSolverTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( DefaultCriterionTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    ContextPtr context   = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );

    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );

    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    DenseVector<ValueType> rhs( coefficients.getColDistributionPtr(), 1.0 );
    rhs.setContextPtr( context );

    DenseVector<ValueType> solution( rhs );
    solution.setContextPtr( context );
    solution.redistribute( coefficients.getRowDistributionPtr() );

    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        Solver* solver( Solver::create( values[i], "" ) );

        IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>( solver );

        if( iterativeSolver == NULL )
        {
            SCAI_LOG_INFO( logger, "Skipping solver " << values[i] << ": no iterative Solver");
            continue;
        }

        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );

        iterativeSolver->initialize( coefficients );
        iterativeSolver->solve( solution, rhs );

        BOOST_CHECK_EQUAL( iterativeSolver->getIterationCount(), 1 );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( IterationCountStoppingCriterionTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "IterativeSolverTest uses context = " << context->getType() );

    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );

    const ValueType solutionInitValue = 1.0;
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), solutionInitValue );
    // TODO: use constructor to set context
    solution.setContextPtr( context );


    DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), solutionInitValue+1.0 );
    // TODO: use constructor to set context
    exactSolution.setContextPtr( context );

    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType numIterations = 300;
    CriterionPtr criterion( new IterationCount( numIterations ) );

    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        Solver* solver = Solver::create( values[i], "" );
        IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>( solver );

        if( iterativeSolver == NULL )
        {
            SCAI_LOG_INFO( logger, "Skipping solver " << values[i] << ": no iterative Solver");
            continue;
        }

        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );

        // TODO: this should be tested for ALL solvers (not only iterative)
        // Solver has to be initialized before solve is called
        //BOOST_CHECK_THROW ( {iterativeSolver->solve( solution, rhs );}, scai::common::Exception );

        iterativeSolver->setStoppingCriterion( criterion );
        iterativeSolver->initialize( coefficients );
        iterativeSolver->solve( solution, rhs );

        BOOST_CHECK( numIterations == iterativeSolver->getIterationCount() );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( SolveTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType N1 = 10;
    const IndexType N2 = 10;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "IterativeSolverTest uses context = " << context->getType() );

    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );

    const ValueType solutionInitValue = 1.0;
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), solutionInitValue );
    // TODO: use constructor to set context
    solution.setContextPtr( context );


    DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), solutionInitValue+1.0 );
    // TODO: use constructor to set context
    exactSolution.setContextPtr( context );

    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType maxExpectedIterations = 3000;
    CriterionPtr criterion( new IterationCount( maxExpectedIterations ) );

    // Get all available solvers
    std::vector<std::string> values;
    Solver::getCreateValues( values );

    const int numSolvers = (int)values.size();

    for(int i=0; i < numSolvers; i++)
    {
        Solver* solver = Solver::create( values[i], "" );
        IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>( solver );

        if( iterativeSolver == NULL )
        {
            SCAI_LOG_INFO( logger, "Skipping solver " << values[i] << ": no iterative Solver");
            continue;
        }

        SCAI_LOG_INFO( logger, "Testing solver " << values[i] );

        iterativeSolver->setStoppingCriterion( criterion );
        iterativeSolver->initialize( coefficients );
        solution = solutionInitValue;
        iterativeSolver->solve( solution, rhs );

        DenseVector<ValueType> diff( solution - exactSolution );

        Scalar s                  = maxNorm( diff );
        ValueType realMaxNorm     = s.getValue<ValueType>();
        ValueType expectedMaxNorm = 1E-4;

        SCAI_LOG_INFO( logger, "maxNorm of diff = " << s << " = ( solution - exactSolution ) = " << realMaxNorm );

        BOOST_CHECK( realMaxNorm < expectedMaxNorm );


        // Test solve with preconditioner
        SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
        iterativeSolver->setPreconditioner( preconditioner );

        // TODO: this should be tested for ALL preconditioners
        // Solver has to be initialized AFTER a preconditioner is set
        //BOOST_CHECK_THROW ( {iterativeSolver->solve( solution, rhs );}, scai::common::Exception );

        iterativeSolver->initialize( coefficients );
        solution = solutionInitValue;
        iterativeSolver->solve( solution, rhs );

        BOOST_CHECK( maxExpectedIterations == iterativeSolver->getIterationCount() );

        diff = solution - exactSolution;

        s               = maxNorm( diff );
        realMaxNorm     = s.getValue<ValueType>();
        expectedMaxNorm = 1E-4;

        SCAI_LOG_INFO( logger, "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << realMaxNorm );

        BOOST_CHECK( realMaxNorm < expectedMaxNorm );
    }
}

BOOST_AUTO_TEST_SUITE_END();

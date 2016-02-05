/**
 * @file P_JacobiTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Contains the implementation of the class P_JacobiTest.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 27.02.2012
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/DefaultJacobi.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr mComm;

struct P_JacobiTestConfig
{
    P_JacobiTestConfig()
    {
        mComm = Communicator::getCommunicator( scai::lama::communicator::MPI );
    }

    ~P_JacobiTestConfig()
    {
        mComm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_JacobiTest, P_JacobiTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_JacobiTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionMethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    DefaultJacobi jacobiSolver( "JacobiTestSolver" );

    if ( SCAI_LOG_INFO_ON( logger ) )
    {
        LoggerPtr slogger( new CommonLogger(
                               "<Jacobi>: ",
                               scai::solver::LogLevel::solverInformation,
                               scai::solver::LoggerWriteBehaviour::toConsoleOnly ) );
        jacobiSolver.setLogger( slogger );
    }

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    // created matrix must have diagonal property
    BOOST_REQUIRE( helpcoefficients.hasDiagonalProperty() );
    MatrixType coefficients( helpcoefficients );
    // converted matrix must have diagonal property
    BOOST_REQUIRE( coefficients.hasDiagonalProperty() );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), mComm ) );
    coefficients.redistribute( dist, dist );
    // converted redistributed matrix must have kept the diagonal property
    BOOST_REQUIRE( coefficients.hasDiagonalProperty() );
    coefficients.setContextPtr( loc );
    DenseVector<ValueType> solution( dist, 1.0 );
    const DenseVector<ValueType> exactSolution( dist, 2.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;
    //initialize
    IndexType expectedIterations = 100;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.initialize( coefficients );
    SCAI_LOG_INFO( logger, "jacobiSolver::coefficients = " << coefficients );
    jacobiSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, jacobiSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "max norm ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithoutPreconditionMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< ELLSparseMatrix<ValueType> >( context );
        // testSolveWithoutPreconditionMethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< COOSparseMatrix<ValueType> >( context );
        // testSolveWithoutPreconditionMethod< DenseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionMethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    DefaultJacobi jacobiSolver( "JacobiTestSolver" );

    if ( SCAI_LOG_INFO_ON( logger ) )
    {
        LoggerPtr slogger( new CommonLogger(
                               "<Jacobi + Precond>: ",
                               scai::solver::LogLevel::solverInformation,
                               scai::solver::LoggerWriteBehaviour::toConsoleOnly ) );
        jacobiSolver.setLogger( slogger );
    }

    CSRSparseMatrix<ValueType> helpcoefficients;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), mComm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContextPtr( loc );
    DenseVector<ValueType> solution( dist, 1.0 );
    const DenseVector<ValueType> exactSolution( dist, 2.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;
    IndexType expectedIterations = 100;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    jacobiSolver.setStoppingCriterion( criterion );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    jacobiSolver.setPreconditioner( preconditioner );
    jacobiSolver.initialize( coefficients );
    jacobiSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, jacobiSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "max norm ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithPreconditionMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< COOSparseMatrix<ValueType> >( context );
        //@todo: DIA has problems with diagonal property after redistribute
        // testSolveWithPreconditionMethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< JDSSparseMatrix<ValueType> >( context );
        //@todo: Dense cannot be constructed by SparseMatrix
        // testSolveWithPreconditionMethod< DenseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

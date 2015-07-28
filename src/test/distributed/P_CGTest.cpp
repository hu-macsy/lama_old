/**
 * @file P_CGTest.cpp
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
 * @brief Contains the implementation of the class P_CGTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 27.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/TrivialPreconditioner.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/TestMacros.hpp>

using namespace lama;
using namespace memory;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_CGTestConfig
{
    P_CGTestConfig()
    {
        comm = Communicator::get( "MPI" );
    }

    ~P_CGTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_CGTest, P_CGTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.P_CGTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    LAMA_LOG_INFO( logger, "testSolveWithoutPreconditionmethod<" << typeid( MatrixType ).name() << " at " << *loc );
    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CG cgSolver( "CGTestSolver" );
    LAMA_LOG_INFO( logger, "Solver = " << cgSolver );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    LAMA_LOG_INFO( logger, "Poisson2D matrix = " << helpcoefficients );
    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "Poisson2D matrix (converted to MatrixType)  = " << helpcoefficients );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContext( loc );
    DenseVector<ValueType> solution( dist, 2.0 );
    const DenseVector<ValueType> exactSolution( dist, 1.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;
    //initialize
    IndexType expectedIterations = 15;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    ValueType sval = s.getValue<ValueType>();

    if ( ! ( sval < 1E-6 ) )
    {
        LAMA_LOG_ERROR( logger, "max norm of diff = " << sval << ", should be < 1E-6 " )
    }

    BOOST_CHECK( sval < 1E-6 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithoutPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< COOSparseMatrix<ValueType> >( context );
        // @todo: does not work as DenseMatrix = SparseMatrix not supported yet
        // testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    LAMA_LOG_INFO( logger, "testSolveWithPreconditionmethod<" << typeid( MatrixType ).name() << "> on " << *loc );
    CG cgSolver( "CGTestSolver" );

    if ( LAMA_LOG_INFO_ON( logger ) )
    {
        LoggerPtr slogger( new CommonLogger(
                               "<SOR>: ",
                               lama::LogLevel::solverInformation,
                               lama::LoggerWriteBehaviour::toConsoleOnly,
                               new Timer() ) );
        cgSolver.setLogger( slogger );
    }

    const IndexType N1 = 4;
    const IndexType N2 = 4;
    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> csrCoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( csrCoefficients, 9, N1, N2 );
    MatrixType coefficients( csrCoefficients );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContext( loc );
    DenseVector<ValueType> solution( dist, 1.0 );
    const DenseVector<ValueType> exactSolution( dist, 2.0 );
    DenseVector<ValueType> rhs( dist, 0.0 );
    rhs = coefficients * exactSolution;
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    cgSolver.setStoppingCriterion( criterion );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    cgSolver.setPreconditioner( preconditioner );
    LAMA_LOG_INFO( logger, "matrix for CG solver = " << coefficients );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "max norm ( solution - exactSolution ) = " << s );

    if ( s.getValue<ValueType>() >= 1E-6 )
    {
        LAMA_LOG_ERROR( logger, "cgSolver for " << coefficients
                        << ": max norm ( solution - exactSolution ) = " << s );
    }

    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< COOSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< DenseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

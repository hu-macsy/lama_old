/**
 * @file P_JacobiTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/DefaultJacobi.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/TrivialPreconditioner.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr mComm;

struct P_JacobiTestConfig
{
    P_JacobiTestConfig()
    {
        mComm = CommunicatorFactory::get( "MPI" );

    }

    ~P_JacobiTestConfig()
    {
        mComm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_JacobiTest, P_JacobiTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.P_JacobiTest" );

/* --------------------------------------------------------------------- */

template<typename mt>
void testSolveWithoutPreconditionMethod( ContextPtr loc )
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    DefaultJacobi jacobiSolver( "JacobiTestSolver" );

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

    coefficients.setContext( loc );

    DenseVector<ValueType> solution( dist, 2.0 );

    const DenseVector<ValueType> exactSolution( dist, 1.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;

    //initialize
    IndexType expectedIterations = 100;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.initialize( coefficients );

    LAMA_LOG_INFO( logger, "jacobiSolver::coefficients = " << coefficients );

    jacobiSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations, jacobiSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "max norm ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, T, test_types ) {
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithoutPreconditionMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< ELLSparseMatrix<ValueType> >( context );
        //testSolveWithoutPreconditionMethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionMethod< COOSparseMatrix<ValueType> >( context );
        //testSolveWithoutPreconditionMethod< DenseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

template<typename mt>
void testSolveWithPreconditionMethod( ContextPtr loc )
{
    typedef typename mt::ValueType ValueType;

//    LoggerPtr slogger( new CommonLogger(
//        "<SOR>: ",
//        lama::LogLevel::solverInformation,
//        lama::LoggerWriteBehaviour::toConsoleOnly,
//        std::auto_ptr<Timer>( new Timer() ) ) );

    DefaultJacobi jacobiSolver( "JacobiTestSolver"/*, slogger */);

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );

    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), mComm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContext( loc );

    DenseVector<ValueType> solution( dist, 1.0 );
    const DenseVector<ValueType> exactSolution( dist, 2.0 );
    DenseVector<ValueType> rhs( dist, 0.0 );
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
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );

}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, T, test_types ) {
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithPreconditionMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< COOSparseMatrix<ValueType> >( context );
        //@todo: DIA has problems with diagonal property after redistribute
        //testSolveWithPreconditionMethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithPreconditionMethod< JDSSparseMatrix<ValueType> >( context );
        //@todo: Dense cannot be constructed by SparseMatrix
        testSolveWithPreconditionMethod< DenseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/**
 * @file P_SORTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the class P_SORTest.
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/SOR.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
//#include <lama/solver/logging/CommonLogger.hpp>
//#include <lama/solver/logging/OpenMPTimer.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/Communicator.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<float,double> test_types;
//typedef boost::mpl::list<float>  test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_SORTestConfig
{
    P_SORTestConfig()
    {
        comm = CommunicatorFactory::get();
    }

    ~P_SORTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( P_SORTest, P_SORTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.P_SORTest" );

/* --------------------------------------------------------------------- */
template<typename mt>
void testSolveMethod( ContextPtr loc )
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    SOR sor( "SORTestSolver" );

    const IndexType N1 = 8;
    const IndexType N2 = 8;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    LAMA_LOG_INFO( logger, "Poisson 2D CSR matrix = " << helpcoefficients );

    MatrixType coefficients( helpcoefficients );

    // redistribution of coefficient: does not keep diagonal property for DIA

    DistributionPtr dist( new BlockDistribution( helpcoefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContext( loc );

    LAMA_LOG_INFO( logger, "SOR: coefficients = " << coefficients );

    const DenseVector<ValueType> exactSolution( dist, static_cast<ValueType>( 1.1 ) );
    DenseVector<ValueType> rhs( dist );
    rhs = coefficients * exactSolution;
    DenseVector<ValueType> solution( dist, 3.0 );

    sor.initialize( coefficients );

    IndexType expectedIterations = 120;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    sor.setStoppingCriterion( criterion );
    sor.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, sor.getIterationCount() );

    DenseVector<ValueType> diff( exactSolution );
    diff = solution - exactSolution;
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm ( solution - exactSolutino ) = " << s );

    //TODO: Check bad value of MaxNorm 0.65
    BOOST_CHECK( s.getValue<ValueType>() < 0.65 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( solveTest, T, test_types ) {
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveMethod< JDSSparseMatrix<ValueType> >( context );
        // @todo: does not work as DIA does not preserve diagonal property
        // testSolveMethod< DIASparseMatrix<ValueType> >( context );
        testSolveMethod< COOSparseMatrix<ValueType> >( context );
        // @todo: not yet available: DenseMatrix = SparseMatrix
        // testSolveMethod< DenseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

template<typename mt>
void testSolvePoissonMethod()
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    double omega;
    double omegaMin = 0.3;
    double omegaMax = 0.7;
    double omegaStep = 0.1;

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    int n = static_cast<int>( ( omegaMax - omegaMin ) / omegaStep + 1 );

    n = 1;

    for ( int i = 0; i < n; i++ )
    {
        omega = omegaMin + omegaStep * i;

        LAMA_LOG_INFO( logger, "run with omega = " << omega );

        //Optional Logger; deactivated for testruns
//        LoggerPtr slogger( new CommonLogger(
//            "<SOR>: ",
//            lama::LogLevel::solverInformation,
//            lama::LoggerWriteBehaviour::toConsoleOnly,
//            std::auto_ptr<Timer>( new OpenMPTimer() ) ) );

        SOR sor( "SORTest", omega /*, slogger*/);

        // create Poisson2D matrix, only for CSR format available

        CSRSparseMatrix<ValueType> csrCoefficients;
        MatrixCreator<ValueType>::buildPoisson2D( csrCoefficients, 9, N1, N2 );

        // copy constructur, converts to format of MatrixType

        MatrixType coefficients( csrCoefficients );

        LAMA_LOG_INFO( logger, "coefficients for SOR: " << csrCoefficients );

        // DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm) );
        // coefficients.redistribute( dist, dist );

        DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 3.0 );
        DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 1.0 );

        DenseVector<ValueType> rhs( coefficients.getDistributionPtr(), 1.0 );

        rhs = coefficients * exactSolution;
        sor.initialize( coefficients );

        //IterationCount
        IndexType expectedIterations = 100;
        CriterionPtr criterion( new IterationCount( expectedIterations ) );
        sor.setStoppingCriterion( criterion ); //what kind of criterion is possible? How can I implement it?

        sor.solve( solution, rhs );

        BOOST_CHECK_EQUAL( expectedIterations, sor.getIterationCount() );
        //Asserts that expectedIterations and made iterations are equals.

        DenseVector<ValueType> diff = solution - exactSolution;
        Scalar solutionPrecision = maxNorm( diff );
        BOOST_CHECK( solutionPrecision.getValue<ValueType>() < 1E-2 );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolvePoisson, T, test_types ) {
// @todo: Test does not work with general distribution
// testSolvePoissonMethod< CSRSparseMatrix<T> >();
// testSolvePoissonMethod< ELLSparseMatrix<T> >();
// testSolvePoissonMethod< JDSSparseMatrix<T> >();
// testSolvePoissonMethod< DIASparseMatrix<T> >();
// testSolvePoissonMethod< COOSparseMatrix<T> >();
}

///* --------------------------------------------------------------------- */

template<typename mt>
void testSolve2Method( ContextPtr loc )
{
    typedef mt MatrixType;
    typedef typename mt::ValueType ValueType;

    double omega;
    double omegaMin = 0.8;
    double omegaMax = 1.0;
    double omegaStep = 0.05;

    int n = static_cast<int>( std::floor( ( omegaMax - omegaMin ) / omegaStep + 1 ) );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    for ( int i = 0; i < n; i++ )
    {
        omega = omegaMin + omegaStep * i;

        LAMA_LOG_INFO( logger, "iterate i = " << i << " of " << n << " iterations" << ", omega = " << omega );

//        LoggerPtr slogger( new CommonLogger(
//            "<SOR>: ",
//            lama::LogLevel::solverInformation,
//            lama::LoggerWriteBehaviour::toConsoleOnly,
//            std::auto_ptr<Timer>( new OpenMPTimer() ) ) );

        SOR sor( "SORTest", omega/*, slogger */);

        CSRSparseMatrix<ValueType> csrCoefficients;
        MatrixCreator<ValueType>::buildPoisson2D( csrCoefficients, 9, N1, N2 );
        DistributionPtr dist( new BlockDistribution( csrCoefficients.getNumRows(), comm ) );
        csrCoefficients.redistribute( dist, dist );

        // Attention: csrCoefficients has a general distribution, so redistribution to
        // block distribution is not really necessary, but more efficient.
        // Nevertheless: hangs with 3 processors if we do not redistribute

        // Attention2: we could apply redistrubte only to coefficients and not to
        // to csrCoefficients. But then the DIA format has problems as it does
        // not keep the diagonal property.

        MatrixType coefficients( csrCoefficients );
        coefficients.setContext( loc );

        DenseVector<ValueType> rhs( dist );
        DenseVector<ValueType> exactSolution( dist, 1.0 );
        DenseVector<ValueType> solution( dist, 2.0 );

        rhs = coefficients * exactSolution;

        sor.initialize( coefficients );

        IndexType expectedIterations = 25;
        CriterionPtr criterion( new IterationCount( expectedIterations ) );
        sor.setStoppingCriterion( criterion );

        //SOLVING
        sor.solve( solution, rhs );
        BOOST_CHECK_EQUAL( expectedIterations, sor.getIterationCount() );

        DenseVector<ValueType> diff = solution - exactSolution;

        Scalar solutionPrecision = maxNorm( diff );

        BOOST_CHECK( solutionPrecision.getValue<ValueType>() < 0.01 );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve2, T, test_types ) {
    typedef T ValueType;

    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolve2Method< CSRSparseMatrix<ValueType> >( context );
        testSolve2Method< ELLSparseMatrix<ValueType> >( context );
        testSolve2Method< JDSSparseMatrix<ValueType> >( context );
        // testSolve2Method< DIASparseMatrix<ValueType> >( context );
        testSolve2Method< COOSparseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

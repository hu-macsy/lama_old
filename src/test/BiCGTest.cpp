/**
 * @file BiCGTest.cpp
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
 * @brief BiCGTest.cpp
 * @author Lauretta Schubert
 * @date 04.07.2013
 * @since 1.1.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/BiCG.hpp>
#include <lama/solver/TrivialPreconditioner.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/logger/Timer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/norm/MaxNorm.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/TestMacros.hpp>

using namespace lama;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BiCGTest );

LAMA_LOG_DEF_LOGGER( logger, "Test.BiCGTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<BiCG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          new Timer() ) );
    BiCG bicgSolver( "BiCGTestSolver", slogger );
    BOOST_CHECK_EQUAL( bicgSolver.getId(), "BiCGTestSolver" );
    BiCG bicgSolver2( "BiCGTestSolver2" );
    BOOST_CHECK_EQUAL( bicgSolver2.getId(), "BiCGTestSolver2" );
    BiCG bicgSolver3( bicgSolver2 );
    BOOST_CHECK_EQUAL( bicgSolver3.getId(), "BiCGTestSolver2" );
    BOOST_CHECK( bicgSolver3.getPreconditioner() == 0 );
    BiCG bicgSolver4( "BiCGTestSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    bicgSolver4.setPreconditioner( preconditioner );
    CriterionPtr criterion( new IterationCount( 10 ) );
    bicgSolver4.setStoppingCriterion( criterion );
    BiCG bicgSolver5( bicgSolver4 );
    BOOST_CHECK_EQUAL( bicgSolver5.getId(), bicgSolver4.getId() );
    BOOST_CHECK_EQUAL( bicgSolver5.getPreconditioner()->getId(), bicgSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    LoggerPtr slogger(
        new CommonLogger( "<BiCG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          new Timer() ) );
    BiCG bicgSolver( "BiCGTestSolver", slogger );
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    // convert to the corresponding matrix type, keep distribution
    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "BiCGTest uses context = " << context->getType() );
    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    bicgSolver.setStoppingCriterion( criterion );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    bicgSolver.setPreconditioner( preconditioner );
    bicgSolver.initialize( coefficients );
    bicgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, bicgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger,
                   "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-4 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< COOSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithPreconditionmethod< DenseMatrix<ValueType> >( context );
        // ToDo: does not work with NP=2:    testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >();
        // ToDo: does not work with NP=2:    testSolveWithPreconditionmethod< DenseMatrix<ValueType> >();
    }
}

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    BiCG bicgSolver( "BiCGTestSolver" );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficient matrix = " << coefficients );
    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "BiCGTest uses context = " << context->getType() );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );
    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()
    const DenseVector<ValueType> rhs( coefficients * exactSolution );
    LAMA_LOG_INFO( logger, "rhs = " << rhs );
    //initialize
    IndexType expectedIterations = 20;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    bicgSolver.setStoppingCriterion( criterion );
    bicgSolver.initialize( coefficients );
    bicgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, bicgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveWithoutPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< COOSparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< DIASparseMatrix<ValueType> >( context );
        testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >( context );
        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DenseMatrix<T> >();
        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DIASparseMatrix<T> >();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE ( simpleTest, ValueType, test_types )
{
    const IndexType n = 3;
    const IndexType numValues = 5;
    IndexType ia[] = { 0, 2, 3, 5 };
    IndexType ja[] = { 0, 2, 1, 1, 2 };
    ValueType matrixValues[] = { 1.0, 2.0, 1.0, 2.0, 1.0 };
    const LAMAArray<IndexType> matrixIA = LAMAArray<IndexType>( n + 1, ia );
    const LAMAArray<IndexType> matrixJA = LAMAArray<IndexType>( numValues, ja );
    const LAMAArray<ValueType> mValues  = LAMAArray<ValueType>( numValues, matrixValues );
    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( n, n, numValues, matrixIA, matrixJA, mValues );
    lama::DistributionPtr dist( new lama::NoDistribution( n ) );
    CSRSparseMatrix<ValueType> matrix( *csrStorage, dist, dist );
    ValueType vectorValues[] = { 3.0f, 0.5f, 2.0f };
    const LAMAArray<ValueType> vValues  = LAMAArray<ValueType>( n, vectorValues );
    DenseVector<ValueType> rhs( n, 0.0 );
    rhs.setValues( vValues );
    ValueType vectorValues2[] = { 1.0f, 0.5f, 1.0f };
    const LAMAArray<ValueType> vValues2  = LAMAArray<ValueType>( n, vectorValues2 );
    DenseVector<ValueType> exactSolution( n, 0.0 );
    exactSolution.setValues( vValues2 );
    DenseVector<ValueType> solution( n, 0.0 );
    //initialize
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    BiCG bicgSolver( "BiCGTestSolver" );
    bicgSolver.setStoppingCriterion( criterion );
    bicgSolver.initialize( matrix );
    bicgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, bicgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE ( testDefaultCriterionSet, ValueType, test_types )
{
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    BiCG bicgSolver( "BiCGTestSolver" );
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 0.0 );
    bicgSolver.initialize( coefficients );
    bicgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( bicgSolver.getIterationCount(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    BiCG bicgSolver( "BiCGTestSolver" );
    LAMA_WRITEAT_TEST( bicgSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    BiCG bicgSolver1( "BiCGTestSolver" );
    SolverPtr solverptr = bicgSolver1.copy();
    BOOST_CHECK_EQUAL( solverptr->getId(), "BiCGTestSolver" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

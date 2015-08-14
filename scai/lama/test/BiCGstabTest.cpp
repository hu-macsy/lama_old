/**
 * @file BiCGstabTest.cpp
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
 * @brief BiCGstabTest.cpp
 * @author lschubert
 * @date 07.08.2013
 * @since 1.1.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/solver/BiCGstab.hpp>
#include <scai/lama/solver/TrivialPreconditioner.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/solver/criteria/ResidualThreshold.hpp>
#include <scai/lama/solver/logger/Timer.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <test/TestMacros.hpp>

using namespace lama;
using namespace memory;

typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BiCGstabTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.BiCGstabTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          new Timer() ) );

    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );
    BOOST_CHECK_EQUAL( BiCGstabSolver.getId(), "BiCGstabTestSolver" );

    BiCGstab BiCGstabSolver2( "BiCGstabTestSolver2" );
    BOOST_CHECK_EQUAL( BiCGstabSolver2.getId(), "BiCGstabTestSolver2" );

    BiCGstab BiCGstabSolver3( BiCGstabSolver2 );
    BOOST_CHECK_EQUAL( BiCGstabSolver3.getId(), "BiCGstabTestSolver2" );
    BOOST_CHECK( BiCGstabSolver3.getPreconditioner() == 0 );

    BiCGstab BiCGstabSolver4( "BiCGstabSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    BiCGstabSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    BiCGstabSolver4.setStoppingCriterion( criterion );

    BiCGstab BiCGstabSolver5( BiCGstabSolver4 );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getId(), BiCGstabSolver4.getId() );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getPreconditioner()->getId(), BiCGstabSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    BiCGstab BiCGstabSolver( "TestBiCGstab" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );

    DenseVector<ValueType> solution( rhs );

    BiCGstabSolver.initialize( coefficients );   // Not WORKING

    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( BiCGstabSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    LoggerPtr slogger(
        new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          new Timer() ) );

    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    // convert to the corresponding matrix type, keep distribution

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficients matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "BiCGstabTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType expectedIterations = 300;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    BiCGstabSolver.setStoppingCriterion( criterion );

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    BiCGstabSolver.setPreconditioner( preconditioner );

    BiCGstabSolver.initialize( coefficients );
    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK( expectedIterations >= BiCGstabSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger,
                   "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-4 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, T, test_types ) {
    typedef T ValueType;

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

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

        LoggerPtr slogger(
        new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          new Timer() ) );

    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );


    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "BiCGstabTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );

    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()

    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    LAMA_LOG_INFO( logger, "rhs = " << rhs );

    //initialize
    IndexType expectedIterations = 300  ;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    BiCGstabSolver.setStoppingCriterion( criterion );
    BiCGstabSolver.initialize( coefficients );

    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK( expectedIterations >= BiCGstabSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-4 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, test_types ) {

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
    ValueType vectorValues[] = { 0.3f, 0.7f, 3.1416f };
    const LAMAArray<ValueType> vValues  = LAMAArray<ValueType>( n, vectorValues );
    DenseVector<ValueType> rhs( n, 0.0 );
    rhs.setValues( vValues );
    ValueType vectorValues2[] = {  -3.183185307179586, 0.7, 1.741592653589793 };
    const LAMAArray<ValueType> vValues2  = LAMAArray<ValueType>( n, vectorValues2 );
    DenseVector<ValueType> exactSolution( n, 0.0 );
    exactSolution.setValues( vValues2 );
    DenseVector<ValueType> solution( n, 0.0 );
    //initialize
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    BiCGstab bicgstabSolver( "BiCGstabTestSolver" );
    bicgstabSolver.setStoppingCriterion( criterion );
    bicgstabSolver.initialize( matrix );
    bicgstabSolver.solve( solution, rhs );
    BOOST_CHECK( expectedIterations >= bicgstabSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-4 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    BiCGstab BiCGstabSolver( "BiCGstabTestSolver" );
    LAMA_WRITEAT_TEST( BiCGstabSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    BiCGstab BiCGstabSolver1( "BiCGTestSolver" );

    SolverPtr solverptr = BiCGstabSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "BiCGTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();











/**
 * @file QMRTest.cpp
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
 * @brief QMRTest.cpp
 * @author lschubert
 * @date 07.08.2013
 * @since 1.1.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/QMR.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

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

#include <scai/lama/test/TestMacros.hpp>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;

typedef boost::mpl::list<double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( QMRTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.QMRTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<QMR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    QMR QMRSolver( "QMRTestSolver", slogger );
    BOOST_CHECK_EQUAL( QMRSolver.getId(), "QMRTestSolver" );

    QMR QMRSolver2( "QMRTestSolver2" );
    BOOST_CHECK_EQUAL( QMRSolver2.getId(), "QMRTestSolver2" );

    QMR QMRSolver3( QMRSolver2 );
    BOOST_CHECK_EQUAL( QMRSolver3.getId(), "QMRTestSolver2" );
    BOOST_CHECK( QMRSolver3.getPreconditioner() == 0 );

    QMR QMRSolver4( "QMRSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    QMRSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    QMRSolver4.setStoppingCriterion( criterion );

    QMR QMRSolver5( QMRSolver4 );
    BOOST_CHECK_EQUAL( QMRSolver5.getId(), QMRSolver4.getId() );
    BOOST_CHECK_EQUAL( QMRSolver5.getPreconditioner()->getId(), QMRSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    QMR QMRSolver( "TestQMR" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );

    DenseVector<ValueType> solution( rhs );

    QMRSolver.initialize( coefficients );   // Not WORKING

    QMRSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( QMRSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    LoggerPtr slogger(
        new CommonLogger( "<QMR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly) );
    QMR QMRSolver( "QMRTestSolver", slogger );
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 5, N1, N2 );
    // convert to the corresponding matrix type, keep distribution
    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "QMRTest uses context = " << context->getType() );
    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );
    IndexType expectedIterations = 5;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    QMRSolver.setStoppingCriterion( criterion );
 
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    preconditioner->initialize(coefficients);
    QMRSolver.setPreconditioner( preconditioner );
    
    QMRSolver.initialize( coefficients );
    QMRSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, QMRSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger,
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
        new CommonLogger( "<QMR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    QMR QMRSolver( "QMRTestSolver", slogger );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );


    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 5, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "QMRTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );
    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    SCAI_LOG_INFO( logger, "rhs = " << rhs );

    IndexType expectedIterations = 5;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    QMRSolver.setStoppingCriterion( criterion );
    QMRSolver.initialize( coefficients );

    QMRSolver.solve( solution, rhs );

    BOOST_CHECK( expectedIterations >= QMRSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
	BOOST_CHECK_SMALL( s.getValue<ValueType>(), common::TypeTraits<ValueType>::small() );
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
    const LArray<IndexType> matrixIA( n + 1, ia );
    const LArray<IndexType> matrixJA( numValues, ja );
    const LArray<ValueType> mValues( numValues, matrixValues );
    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( n, n, numValues, matrixIA, matrixJA, mValues );
    scai::lama::DistributionPtr dist( new scai::lama::NoDistribution( n ) );
    CSRSparseMatrix<ValueType> matrix( *csrStorage, dist, dist );
    ValueType vectorValues[] = { 0.3f, 0.7f, 3.1416f };
    const LArray<ValueType> vValues ( n, vectorValues );
    DenseVector<ValueType> rhs( n, 0.0 );
    rhs.setValues( vValues );
    ValueType vectorValues2[] = {  -3.183185307179586, 0.7, 1.741592653589793 };
    const LArray<ValueType> vValues2( n, vectorValues2 );
    DenseVector<ValueType> exactSolution( n, 0.0 );
    exactSolution.setValues( vValues2 );
    DenseVector<ValueType> solution( n, 0.0 );
    //initialize
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    QMR QMRSolver( "QMRTestSolver" );
    QMRSolver.setStoppingCriterion( criterion );
    QMRSolver.initialize( matrix );
    QMRSolver.solve( solution, rhs );
    BOOST_CHECK( expectedIterations >= QMRSolver.getIterationCount() );

	SCAI_LOG_INFO( logger, "number of iterations " << QMRSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
   // BOOST_CHECK( s.getValue<ValueType>() < eps<ValueType>() );
//	BOOST_CHECK_SMALL( s.getValue<ValueType>(), eps<ValueType>() );
	BOOST_CHECK_SMALL( s.getValue<ValueType>(), ValueType(1e-4) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    QMR QMRSolver( "QMRTestSolver" );
    LAMA_WRITEAT_TEST( QMRSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    QMR QMRSolver1( "QMRTestSolver" );

    SolverPtr solverptr = QMRSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "QMRTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();











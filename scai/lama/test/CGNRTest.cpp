#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/solver/CGNR.hpp>
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

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CGNRTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGNRTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<CGNR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGNR CGNRSolver( "CGNRTestSolver", slogger );
    BOOST_CHECK_EQUAL( CGNRSolver.getId(), "CGNRTestSolver" );

    CGNR CGNRSolver2( "CGNRTestSolver2" );
    BOOST_CHECK_EQUAL( CGNRSolver2.getId(), "CGNRTestSolver2" );

    CGNR CGNRSolver3( CGNRSolver2 );
    BOOST_CHECK_EQUAL( CGNRSolver3.getId(), "CGNRTestSolver2" );
    BOOST_CHECK( CGNRSolver3.getPreconditioner() == 0 );

    CGNR CGNRSolver4( "CGNRSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGNRSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    CGNRSolver4.setStoppingCriterion( criterion );

    CGNR CGNRSolver5( CGNRSolver4 );
    BOOST_CHECK_EQUAL( CGNRSolver5.getId(), CGNRSolver4.getId() );
    BOOST_CHECK_EQUAL( CGNRSolver5.getPreconditioner()->getId(), CGNRSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    CGNR CGNRSolver( "TestCGNR" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );

    DenseVector<ValueType> solution( rhs );

    CGNRSolver.initialize( coefficients );   // Not WORKING

    CGNRSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( CGNRSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    LoggerPtr slogger(
        new CommonLogger( "<CGNR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGNR CGNRSolver( "CGNRTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    // convert to the corresponding matrix type, keep distribution

    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );

    coefficients.setContext( context );
    SCAI_LOG_INFO( logger, "CGNRTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType expectedIterations = 200;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    CGNRSolver.setStoppingCriterion( criterion );

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGNRSolver.setPreconditioner( preconditioner );

    CGNRSolver.initialize( coefficients );
    CGNRSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations,CGNRSolver.getIterationCount() );

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
   //     testSolveWithPreconditionmethod< DenseMatrix<ValueType> >( context );

        // ToDo: does not work with NP=2:    testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >();
        // ToDo: does not work with NP=2:    testSolveWithPreconditionmethod< DenseMatrix<ValueType> >();
    }
} 

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

        LoggerPtr slogger(
        new CommonLogger( "<CGNR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGNR CGNRSolver( "CGNRTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );


    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContext( context );
    SCAI_LOG_INFO( logger, "CGNRTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );

    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()

    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    SCAI_LOG_INFO( logger, "rhs = " << rhs );

    //initialize
    IndexType expectedIterations = 200  ;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    CGNRSolver.setStoppingCriterion( criterion );
    CGNRSolver.initialize( coefficients );

    CGNRSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations,CGNRSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
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
    //    testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >( context );

        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DenseMatrix<T> >();
        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DIASparseMatrix<T> >();
    }
}

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    CGNR CGNRSolver( "CGNRTestSolver" );
    LAMA_WRITEAT_TEST( CGNRSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    CGNR CGNRSolver1( "CGNRTestSolver" );

    SolverPtr solverptr = CGNRSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "CGNRTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();



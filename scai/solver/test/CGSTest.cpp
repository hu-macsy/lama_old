#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/CGS.hpp>
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

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;



typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CGSTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGSTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<CGS>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGS CGSSolver( "CGSTestSolver", slogger );
    BOOST_CHECK_EQUAL( CGSSolver.getId(), "CGSTestSolver" );

    CGS CGSSolver2( "CGSTestSolver2" );
    BOOST_CHECK_EQUAL( CGSSolver2.getId(), "CGSTestSolver2" );

    CGS CGSSolver3( CGSSolver2 );
    BOOST_CHECK_EQUAL( CGSSolver3.getId(), "CGSTestSolver2" );
    BOOST_CHECK( CGSSolver3.getPreconditioner() == 0 );

    CGS CGSSolver4( "CGSSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGSSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    CGSSolver4.setStoppingCriterion( criterion );

    CGS CGSSolver5( CGSSolver4 );
    BOOST_CHECK_EQUAL( CGSSolver5.getId(), CGSSolver4.getId() );
    BOOST_CHECK_EQUAL( CGSSolver5.getPreconditioner()->getId(), CGSSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    CGS CGSSolver( "TestCGS" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );

    DenseVector<ValueType> solution( rhs );

    CGSSolver.initialize( coefficients );   // Not WORKING

    CGSSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( CGSSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    LoggerPtr slogger(
        new CommonLogger( "<CGS>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGS CGSSolver( "CGSTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
                
    // convert to the corresponding matrix type, keep distribution

    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );

    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "CGSTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType expectedIterations = 300;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    CGSSolver.setStoppingCriterion( criterion );

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGSSolver.setPreconditioner( preconditioner );

    CGSSolver.initialize( coefficients );
    CGSSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations, CGSSolver.getIterationCount() );

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
        new CommonLogger( "<CGS>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGS CGSSolver( "CGSTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );


    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "CGSTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );

    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()

    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    SCAI_LOG_INFO( logger, "rhs = " << rhs );

    //initialize
    IndexType expectedIterations = 300  ;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    CGSSolver.setStoppingCriterion( criterion );
    CGSSolver.initialize( coefficients );

    CGSSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations ,CGSSolver.getIterationCount() );

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
        testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >( context );
        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DenseMatrix<T> >();
        // ToDo: does not run for NP=2: testSolveWithoutPreconditionmethod< DIASparseMatrix<T> >();
    }
}

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    CGS CGSSolver( "CGSTestSolver" );
    LAMA_WRITEAT_TEST( CGSSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    CGS CGSSolver1( "CGTestSolver" );

    SolverPtr solverptr = CGSSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "CGTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();



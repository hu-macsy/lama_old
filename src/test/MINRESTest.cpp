#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/MINRES.hpp>
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



typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MINRESTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.MINRESTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<MINRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          common::unique_ptr<Timer>( new Timer() ) ) );

    MINRES MINRESSolver( "MINRESTestSolver", slogger );
    BOOST_CHECK_EQUAL( MINRESSolver.getId(), "MINRESTestSolver" );

    MINRES MINRESSolver2( "MINRESTestSolver2" );
    BOOST_CHECK_EQUAL( MINRESSolver2.getId(), "MINRESTestSolver2" );

    MINRES MINRESSolver3( MINRESSolver2 );
    BOOST_CHECK_EQUAL( MINRESSolver3.getId(), "MINRESTestSolver2" );
    BOOST_CHECK( MINRESSolver3.getPreconditioner() == 0 );

    MINRES MINRESSolver4( "MINRESSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    MINRESSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    MINRESSolver4.setStoppingCriterion( criterion );

    MINRES MINRESSolver5( MINRESSolver4 );
    BOOST_CHECK_EQUAL( MINRESSolver5.getId(), MINRESSolver4.getId() );
    BOOST_CHECK_EQUAL( MINRESSolver5.getPreconditioner()->getId(), MINRESSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    MINRES MINRESSolver( "TestMINRES" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );

    DenseVector<ValueType> solution( rhs );

    MINRESSolver.initialize( coefficients );   // Not WORKING

    MINRESSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( MINRESSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    LoggerPtr slogger(
        new CommonLogger( "<MINRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          common::unique_ptr<Timer>( new Timer() ) ) );

    MINRES MINRESSolver( "MINRESTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    // convert to the corresponding matrix type, keep distribution

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficients matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "MINRESTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType expectedIterations = 100;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    MINRESSolver.setStoppingCriterion( criterion );

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    MINRESSolver.setPreconditioner( preconditioner );

    MINRESSolver.initialize( coefficients );
    MINRESSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations,MINRESSolver.getIterationCount() );

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
        new CommonLogger( "<MINRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          common::unique_ptr<Timer>( new Timer() ) ) );

    MINRES MINRESSolver( "MINRESTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );


    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "MINRESTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );

    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()

    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    LAMA_LOG_INFO( logger, "rhs = " << rhs );

    //initialize
    IndexType expectedIterations = 100  ;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    MINRESSolver.setStoppingCriterion( criterion );
    MINRESSolver.initialize( coefficients );

    MINRESSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations,MINRESSolver.getIterationCount() );

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

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    MINRES MINRESSolver( "MINRESTestSolver" );
    LAMA_WRITEAT_TEST( MINRESSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    MINRES MINRESSolver1( "MINRESTestSolver" );

    SolverPtr solverptr = MINRESSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "MINRESTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();



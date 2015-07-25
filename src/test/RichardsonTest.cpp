/**
 * @file RichardsonTest.cpp
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
 * @brief Contains the implementation of the class RichardsonTest.
 * @author: 
 * @date 17.04.2015
 * @since 
 **/
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/Richardson.hpp>
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

typedef boost::mpl::list<float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( RichardsonTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.RichardsonTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<Richardson>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );

    Richardson rSolver( "RichardsonTestSolver", slogger );
    BOOST_CHECK_EQUAL( rSolver.getId(), "RichardsonTestSolver" );

    Richardson rSolver2( "RichardsonTestSolver2" );
    BOOST_CHECK_EQUAL( rSolver2.getId(), "RichardsonTestSolver2" );

    Richardson rSolver3( rSolver2 );
    BOOST_CHECK_EQUAL( rSolver3.getId(), "RichardsonTestSolver2" );
    BOOST_CHECK( rSolver3.getPreconditioner() == 0 );

    Richardson rSolver4( "RichardsonTestSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    rSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    rSolver4.setStoppingCriterion( criterion );

    Richardson rSolver5( rSolver4 );
    BOOST_CHECK_EQUAL( rSolver5.getId(), rSolver4.getId() );
    BOOST_CHECK_EQUAL( rSolver5.getPreconditioner()->getId(), rSolver4.getPreconditioner()->getId() );


}
/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;


    LoggerPtr slogger(
        new CommonLogger( "<Richardson>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );

    

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );


    Richardson rSolver( "RichardsonTestSolver" ,slogger ); 
    // convert to the corresponding matrix type, keep distribution

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficients matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "RichardsonTest uses context = " << context->getType() );

    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType expectedIterations = 200;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    rSolver.setStoppingCriterion( criterion );

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    rSolver.setPreconditioner( preconditioner );

    rSolver.initialize( coefficients );
    rSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations, rSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );

    //CAREFUL: Scalar overflow from inf to -nan casted to 0.

    LAMA_LOG_INFO( logger,
                   "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << s.getValue<ValueType>() );

    BOOST_CHECK( s.getValue<ValueType>() < 1E-04 );


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

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    LoggerPtr slogger(
        new CommonLogger( "<Richardson>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );



    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

    MatrixType coefficients( helpcoefficients );
    LAMA_LOG_INFO( logger, "coefficient matrix = " << coefficients );

    coefficients.setContext( context );
    LAMA_LOG_INFO( logger, "RichardsonTest uses context = " << context->getType() );

    Richardson rSolver( "RichardsonTestSolver",slogger);    

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );

    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()

    const DenseVector<ValueType> rhs( coefficients * exactSolution );

    LAMA_LOG_INFO( logger, "rhs = " << rhs );

    //initialize
    IndexType expectedIterations = 200;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    rSolver.setStoppingCriterion( criterion );
    rSolver.initialize( coefficients );

    rSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( expectedIterations, rSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    LAMA_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
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

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    typedef double ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    Richardson rSolver( "RichardsonTestSolver" );

    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );

    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 0.0 );

    rSolver.initialize( coefficients );

    rSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( rSolver.getIterationCount(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Richardson rSolver( "RichardsonTestSolver" );
    LAMA_WRITEAT_TEST( rSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    Richardson rSolver1( "RichardsonTestSolver" );

    SolverPtr solverptr = rSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "RichardsonTestSolver" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

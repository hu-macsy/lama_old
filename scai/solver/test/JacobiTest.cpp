/**
 * @file Jacobi.cpp
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
 * @brief Contains the implementation of the class Jacobi.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 27.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/Jacobi.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

typedef boost::mpl::list<float, double> test_types;

// ---------------------------------------------------------------------------------------------------------------

struct JacobiTestConfig
{
    JacobiTestConfig()
    {
        LoggerPtr loggerD(
            new CommonLogger( "<Jacobi>: ", scai::solver::LogLevel::completeInformation,
                              scai::solver::LoggerWriteBehaviour::toConsoleOnly ) );
        mJacobiDouble = new Jacobi( "JacobiTest double solver", loggerD );
        mJacobiFloat = new Jacobi( "JacobiTest float solver", loggerD );
    }

    ~JacobiTestConfig()
    {
        delete mJacobiDouble;
        delete mJacobiFloat;
    }

    OmegaSolver* mJacobiDouble;
    OmegaSolver* mJacobiFloat;
};

BOOST_FIXTURE_TEST_SUITE( JacobiTest, JacobiTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.JacobiTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testGetId )
{
    BOOST_CHECK_EQUAL( 0, ( *mJacobiDouble ).getId().compare( "JacobiTest double solver" ) );
    BOOST_CHECK_EQUAL( 0, ( *mJacobiFloat ).getId().compare( "JacobiTest float solver" ) );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<SJ>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    Jacobi sjSolver( "SJTestSolver", slogger );
    BOOST_CHECK_EQUAL( sjSolver.getId(), "SJTestSolver" );

    Jacobi sjSolver2( "SJTestSolver2" );
    BOOST_CHECK_EQUAL( sjSolver2.getId(), "SJTestSolver2" );

    Jacobi sjSolver3( sjSolver2 );
    BOOST_CHECK_EQUAL( sjSolver3.getId(), "SJTestSolver2" );
    BOOST_CHECK( sjSolver3.getPreconditioner() == 0 );

    Jacobi sjSolver4( "sjSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    sjSolver4.setPreconditioner( preconditioner );

    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    sjSolver4.setStoppingCriterion( criterion );

    Jacobi sjSolver5( sjSolver4 );
    BOOST_CHECK_EQUAL( sjSolver5.getId(), sjSolver4.getId() );
    BOOST_CHECK_EQUAL( sjSolver5.getPreconditioner()->getId(), sjSolver4.getPreconditioner()->getId() );

    Jacobi sjSolver6( "sjSolver6", 0.05, slogger );
    BOOST_CHECK_EQUAL( sjSolver6.getId(), "sjSolver6" );
    BOOST_CHECK_EQUAL( sjSolver6.getOmega(), 0.05 );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSetAndGetOmega )
{
    ( *mJacobiDouble ).setOmega( 0.6 );
    ( *mJacobiFloat ).setOmega( 0.8 );
    BOOST_CHECK_EQUAL( 0.6, ( *mJacobiDouble ).getOmega().getValue<double>() );
    BOOST_CHECK_EQUAL( 0.8f, ( *mJacobiFloat ).getOmega().getValue<float>() );
}

// ---------------------------------------------------------------------------------------------------------------

template<typename MatrixType>
void testSolveMethod( std::string solverId, ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::string id = solverId;
    LoggerPtr slogger(
        new CommonLogger( solverId, LogLevel::solverInformation, LoggerWriteBehaviour::toConsoleOnly ) );
    Jacobi jacobiSolver( "JacobiTest"/*, slogger*/ );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    CSRSparseMatrix<ValueType> matrix( system.coefficients );

    MatrixType coefficients( matrix );
    coefficients.setContextPtr( context );

    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );

    SCAI_LOG_INFO( logger, "JacobiTest uses context = " << context->getType() );

    DenseVector<ValueType> rhs( system.rhs );
    rhs.setContextPtr( context );
    rhs.redistribute( coefficients.getRowDistributionPtr() );

    DenseVector<ValueType> solution( system.coefficients.getNumRows(), static_cast<ValueType>( 2.1 ) );
    solution.setContextPtr( context );
    solution.redistribute( coefficients.getColDistributionPtr() );

    DenseVector<ValueType> exactSolution( system.solution );
    exactSolution.setContextPtr( context );
    exactSolution.redistribute( coefficients.getColDistributionPtr() );

    jacobiSolver.initialize( coefficients );
    CriterionPtr criterion( new IterationCount( 120 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );

    DenseVector<ValueType> diff( solution - exactSolution );
    diff.setContextPtr( context );
    diff.redistribute( coefficients.getColDistributionPtr() );

    L2Norm l2Norm;
    Scalar norm = l2Norm( diff );
    BOOST_CHECK( norm.getValue<ValueType>() < 1e-1 );
    //bad omega
    //test for even iterations

    DenseVector<ValueType> solutionA( system.coefficients.getNumRows(), 1.0 );
    solutionA.setContextPtr( context );
    solutionA.redistribute( coefficients.getColDistributionPtr() );

    Jacobi jacobiSolverA( "JacobiTest solver 2" );
    jacobiSolverA.initialize( coefficients );
    jacobiSolverA.setOmega( 0.5 );
    jacobiSolverA.solve( solutionA, rhs );
    jacobiSolverA.solve( solutionA, rhs ); //twice

    DenseVector<ValueType> solutionB( system.coefficients.getNumRows(), 1.0 );
    solutionB.setContextPtr( context );
    solutionB.redistribute( coefficients.getColDistributionPtr() );

    Jacobi jacobiSolverB( "JacobiTest solver 2" );
    CriterionPtr criterionB( new IterationCount( 2 ) );
    jacobiSolverB.setStoppingCriterion( criterionB );
    jacobiSolverB.initialize( coefficients );
    jacobiSolverB.setOmega( 0.5 );
    jacobiSolverB.solve( solutionB, rhs );

    DenseVector<ValueType> diffAB( solutionA - solutionB );
    diffAB.setContextPtr( context );
    diffAB.redistribute( coefficients.getColDistributionPtr() );

    Scalar l2norm = l2Norm( diffAB );
    BOOST_CHECK( l2norm.getValue<ValueType>() < 1e-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, test_types )
{
    ContextPtr context = Context::getContextPtr();
    testSolveMethod<CSRSparseMatrix<ValueType> >( "<JacobiCSR> ", context );
    testSolveMethod<ELLSparseMatrix<ValueType> >( "<JacobiELL> ", context );
    testSolveMethod<JDSSparseMatrix<ValueType> >( "<JacobiJDS>", context );
    testSolveMethod<DIASparseMatrix<ValueType> >( "<JacobiDIA>", context );
    testSolveMethod<COOSparseMatrix<ValueType> >( "<JacobiCOO> ", context );
    // TODO: Not working with Dense!
    //testSolveMethod<DenseMatrix<ValueType> >( "<JacobiDense>", context );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
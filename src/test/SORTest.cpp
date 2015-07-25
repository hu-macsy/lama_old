/**
 * @file SORTest.cpp
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
 * @brief Contains the implementation of the class SORTest.
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/solver/SOR.hpp>
#include <lama/solver/TrivialPreconditioner.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/Timer.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>

#include <lama/norm/L2Norm.hpp>
#include <lama/norm/MaxNorm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SORTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.SORTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<SOR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );
    SOR sorSolver( "SORTestSolver", slogger );
    BOOST_CHECK_EQUAL( sorSolver.getId(), "SORTestSolver" );
    SOR sorSolver2( "SORTestSolver2" );
    BOOST_CHECK_EQUAL( sorSolver2.getId(), "SORTestSolver2" );
    SOR sorSolver3( sorSolver2 );
    BOOST_CHECK_EQUAL( sorSolver3.getId(), "SORTestSolver2" );
    BOOST_CHECK( sorSolver3.getPreconditioner() == 0 );
    SOR sorSolver4( "SORTestSolver4", 0.05, slogger );
    BOOST_CHECK_EQUAL( sorSolver4.getOmega(), 0.05 );
    BOOST_CHECK_EQUAL( sorSolver4.getId(), "SORTestSolver4" );
    SOR sorSolver5( "sorSolver5" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner", slogger ) );
    sorSolver5.setPreconditioner( preconditioner );
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    sorSolver5.setStoppingCriterion( criterion );
    SOR sorSolver6( sorSolver5 );
    BOOST_CHECK_EQUAL( sorSolver6.getId(), sorSolver5.getId() );
    BOOST_CHECK_EQUAL( sorSolver6.getOmega(), sorSolver5.getOmega() );
    BOOST_CHECK_EQUAL( sorSolver6.getPreconditioner()->getId(), sorSolver5.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SOROmegaSetandGetTest )
{
    SOR sor( "SOROmegaTest" );
    Scalar omega = 0.05;
    sor.setOmega( omega );
    BOOST_CHECK_EQUAL( sor.getOmega().getValue<double>(), omega.getValue<double>() );
}

/* --------------------------------------------------------------------- */

//testsolve has stoppingCriterion ResidualThreshold
template<typename MatrixType>
void testSolveOmegaMethod( Scalar omega )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    Timer* timer = new Timer();
    std::auto_ptr<Timer> autoTimer( timer );
    LoggerPtr loggerD(
        new CommonLogger( "<SOR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );
    SOR sor( "SORTest", omega, loggerD );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();
    IndexType n = system.rhs.size();
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( n, -1.0 );
    DenseVector<ValueType> correctSolution( system.solution );
    sor.initialize( coefficients );
    NormPtr l2Norm( new L2Norm() );
    CriterionPtr criterionRes( new ResidualThreshold( l2Norm, Scalar( 1e-6 ), ResidualThreshold::Absolute ) );
    sor.setStoppingCriterion( criterionRes );
    sor.solve( solution, rhs );
    DenseVector<ValueType> diff( n, 1.0 );
    diff = solution - correctSolution;
    Scalar solutionPrecision = maxNorm( diff );
    BOOST_CHECK( solutionPrecision.getValue<ValueType>() < 1E-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveOmega, ValueType, test_types )
{
    testSolveOmegaMethod< CSRSparseMatrix<ValueType> >( 0.20 );
    testSolveOmegaMethod< ELLSparseMatrix<ValueType> >( 0.20 );
    testSolveOmegaMethod< JDSSparseMatrix<ValueType> >( 0.20 );
    testSolveOmegaMethod< DIASparseMatrix<ValueType> >( 0.20 );
    testSolveOmegaMethod< COOSparseMatrix<ValueType> >( 0.20 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    SOR sor( "SORTest", 0.05 );
    LAMA_WRITEAT_TEST( sor );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    SOR sorSolver1( "SORTestSolver" );
    SolverPtr solverptr = sorSolver1.copy();
    BOOST_CHECK_EQUAL( solverptr->getId(), "SORTestSolver" );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

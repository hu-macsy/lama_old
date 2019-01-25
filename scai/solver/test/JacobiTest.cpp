/**
 * @file solver/test/JacobiTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Specific tests for the solver class Jacobi.
 * @author Matthias Makulla
 * @date 27.02.2012
 */

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

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;
using namespace dmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( JacobiTest )

// ---------------------------------------------------------------------------------------------------------------

SCAI_LOG_DEF_LOGGER( logger, "Test.JacobiTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( testGetId, ValueType, scai_numeric_test_types )
{
    Jacobi<ValueType> solver( "JacobiTest solver" );
    BOOST_CHECK_EQUAL( 0, solver.getId().compare( "JacobiTest solver" ) );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<SJ>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    Jacobi<ValueType> sjSolver( "SJTestSolver", slogger );
    BOOST_CHECK_EQUAL( sjSolver.getId(), "SJTestSolver" );
    Jacobi<ValueType> sjSolver2( "SJTestSolver2" );
    BOOST_CHECK_EQUAL( sjSolver2.getId(), "SJTestSolver2" );
    Jacobi<ValueType> sjSolver3( sjSolver2 );
    BOOST_CHECK_EQUAL( sjSolver3.getId(), "SJTestSolver2" );
    BOOST_CHECK( sjSolver3.getPreconditioner() == 0 );
    Jacobi<ValueType> sjSolver4( "sjSolver4" );
    SolverPtr<ValueType> preconditioner( new TrivialPreconditioner<ValueType>( "Trivial preconditioner" ) );
    sjSolver4.setPreconditioner( preconditioner );
    IndexType expectedIterations = 10;
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( expectedIterations ) );
    sjSolver4.setStoppingCriterion( criterion );
    Jacobi<ValueType> sjSolver5( sjSolver4 );
    BOOST_CHECK_EQUAL( sjSolver5.getId(), sjSolver4.getId() );
    BOOST_CHECK_EQUAL( sjSolver5.getPreconditioner()->getId(), sjSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSetAndGetOmega )
{
    typedef SCAI_TEST_TYPE ValueType;  // should work fine for all types if it works for one type

    ValueType omega  = 0.05;
    ValueType omega2 = 0.3;

    Jacobi<ValueType> solver( "omegaSolver", omega );
    BOOST_CHECK_EQUAL( solver.getId(), "omegaSolver" );
    BOOST_CHECK_EQUAL( solver.getOmega(), omega );
    solver.setOmega( omega2 );
    BOOST_CHECK_EQUAL( solver.getOmega(), omega2 );
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
    Jacobi<ValueType> jacobiSolver( "JacobiTest"/*, slogger*/ );

    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    CSRSparseMatrix<ValueType> matrix( system.coefficients );
    MatrixType coefficients( matrix );
    coefficients.setContextPtr( context );
    auto dist = std::make_shared<BlockDistribution>( coefficients.getNumRows(), comm );
    coefficients.redistribute( dist, dist );
    SCAI_LOG_INFO( logger, "JacobiTest uses context = " << context->getType() );
    DenseVector<ValueType> rhs( system.rhs );
    rhs.setContextPtr( context );
    rhs.redistribute( coefficients.getRowDistributionPtr() );
    auto solution = denseVectorFill<ValueType>( coefficients.getColDistributionPtr(), ValueType( 2.1 ), context );
    DenseVector<ValueType> exactSolution( system.solution );
    exactSolution.setContextPtr( context );
    exactSolution.redistribute( coefficients.getColDistributionPtr() );
    jacobiSolver.initialize( coefficients );

    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 120 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );
    auto diff = denseVectorEval( solution - exactSolution, context );
    diff.redistribute( coefficients.getColDistributionPtr() );
    L2Norm<ValueType> l2Norm;
    RealType<ValueType> norm = l2Norm( diff );
    RealType<ValueType> eps1 = 1e-1;
    BOOST_CHECK( norm < eps1 );
    //bad omega
    //test for even iterations

    ValueType omega = 0.5;

    auto colDist   = coefficients.getColDistributionPtr();
    auto solutionA = denseVectorFill<ValueType>( colDist, 1.0, context );

    Jacobi<ValueType> jacobiSolverA( "JacobiTest solver 2" );
    jacobiSolverA.initialize( coefficients );
    jacobiSolverA.setOmega( omega );
    jacobiSolverA.solve( solutionA, rhs );
    jacobiSolverA.solve( solutionA, rhs ); //twice

    auto solutionB = denseVectorFill<ValueType>( colDist, 1.0, context );

    Jacobi<ValueType> jacobiSolverB( "JacobiTest solver 2" );
    CriterionPtr<ValueType> criterionB( new IterationCount<ValueType>( 2 ) );
    jacobiSolverB.setStoppingCriterion( criterionB );
    jacobiSolverB.initialize( coefficients );
    jacobiSolverB.setOmega( omega );
    jacobiSolverB.solve( solutionB, rhs );

    auto diffAB = denseVectorEval( solutionA - solutionB, context );

    RealType<ValueType> l2norm = l2Norm( diffAB );
    RealType<ValueType> eps    = 1e-5;
    BOOST_CHECK( l2norm < eps );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, scai_numeric_test_types )
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

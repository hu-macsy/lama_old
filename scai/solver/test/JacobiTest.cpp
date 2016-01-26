/**
 * @file JacobiTest.cpp
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
 * @brief Contains the implementation of the class JacobiTest.cpp
 * @author: Alexander Büchel, Matthias Makulla
 * @date 22.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/DefaultJacobi.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

struct JacobiTestConfig
{
    JacobiTestConfig()
    {
        LoggerPtr loggerD(
            new CommonLogger( "<Jacobi>: ", LogLevel::completeInformation,
                              LoggerWriteBehaviour::toConsoleOnly ) );
        mJacobiDouble = new DefaultJacobi( "JacobiTest double solver", loggerD );
        mJacobiFloat = new DefaultJacobi( "JacobiTest float solver", loggerD );
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

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testDefaultCriterionSet, ValueType, test_types )
{
    DefaultJacobi jacobi( "TestJacobi" );
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    const DenseVector<ValueType> rhs( coefficients.getLocalNumRows(), 1.0 );
    DenseVector<ValueType> solution( rhs );
    jacobi.initialize( coefficients );
    jacobi.solve( solution, rhs );
    BOOST_CHECK_EQUAL( jacobi.getIterationCount(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testSetAndGetOmega )
{
    ( *mJacobiDouble ).setOmega( 0.6 );
    ( *mJacobiFloat ).setOmega( 0.8 );
    BOOST_CHECK_EQUAL( 0.6, ( *mJacobiDouble ).getOmega().getValue<double>() );
    BOOST_CHECK_EQUAL( 0.8f, ( *mJacobiFloat ).getOmega().getValue<float>() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testGetId )
{
    BOOST_CHECK_EQUAL( 0, ( *mJacobiDouble ).getId().compare( "JacobiTest double solver" ) );
    BOOST_CHECK_EQUAL( 0, ( *mJacobiFloat ).getId().compare( "JacobiTest float solver" ) );
}

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveMethod( std::string solverId, ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    std::string id = solverId;
    LoggerPtr slogger(
        new CommonLogger( solverId, LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly) );
    DefaultJacobi jacobiSolver( "JacobiTest", slogger );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    CSRSparseMatrix<ValueType> matrix( system.coefficients );
    MatrixType coefficients( matrix );
    SCAI_LOG_INFO( logger, "coefficient matrix = " << coefficients );
    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "JacobiTest uses context = " << context->getType() );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( system.coefficients.getNumRows(), static_cast<ValueType>( 2.1 ) );
    DenseVector<ValueType> exactSolution( system.solution );
    jacobiSolver.initialize( system.coefficients );
    CriterionPtr criterion( new IterationCount( 120 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );
    DenseVector<ValueType> diff( solution - exactSolution );
    L2Norm l2Norm;
    Scalar norm = l2Norm( diff );
    BOOST_CHECK( norm.getValue<ValueType>() < 1e-1 );
    //bad omega
    //test for even iterations
    DenseVector<ValueType> solutionA( system.coefficients.getNumRows(), 1.0 );
    DefaultJacobi jacobiSolverA( "JacobiTest solver 2" );
    jacobiSolverA.initialize( coefficients );
    jacobiSolverA.setOmega( 0.5 );
    jacobiSolverA.solve( solutionA, rhs );
    jacobiSolverA.solve( solutionA, rhs ); //twice
    DenseVector<ValueType> solutionB( system.coefficients.getNumRows(), 1.0 );
    DefaultJacobi jacobiSolverB( "JacobiTest solver 2" );
    CriterionPtr criterionB( new IterationCount( 2 ) );
    jacobiSolverB.setStoppingCriterion( criterionB );
    jacobiSolverB.initialize( coefficients );
    jacobiSolverB.setOmega( 0.5 );
    jacobiSolverB.solve( solutionB, rhs );
    DenseVector<ValueType> diffAB( solutionA - solutionB );
    Scalar l2norm = l2Norm( diffAB );
    BOOST_CHECK( l2norm.getValue<ValueType>() < 1e-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveMethod<CSRSparseMatrix<ValueType> >( "<JacobiCSR> ", context );
        testSolveMethod<ELLSparseMatrix<ValueType> >( "<JacobiELL> ", context );
        testSolveMethod<JDSSparseMatrix<ValueType> >( "<JacobiJDS>", context );
        testSolveMethod<DIASparseMatrix<ValueType> >( "<JacobiDIA>", context );
        testSolveMethod<COOSparseMatrix<ValueType> >( "<JacobiCOO> ", context );
        testSolveMethod<DenseMatrix<ValueType> >( "<JacobiDense>", context );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    DefaultJacobi jacobiSolver( "JacobiTest" );
    LAMA_WRITEAT_TEST( jacobiSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    DefaultJacobi jacobiSolver1( "JacobiTestSolver" );
    SolverPtr solverptr = jacobiSolver1.copy();
    BOOST_CHECK_EQUAL( solverptr->getId(), "JacobiTestSolver" );
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

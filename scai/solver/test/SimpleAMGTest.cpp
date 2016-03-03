/**
 * @file SimpleAMGTest.cpp
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
 * @brief Contains the implementation of the class SimpleAMGTest.cpp
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SimpleAMGTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SimpleAMGTest" )

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void solverTestMethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    LoggerPtr consoleLogger(
        new CommonLogger( "<SimpleAMG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();
    DenseVector<ValueType> solution( system.coefficients.getNumRows(), 0.0 );
    const DenseVector<ValueType> rhs( system.rhs );
    const DenseVector<ValueType> refSolution( system.solution );
    MatrixType coefficients( system.coefficients );
    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "SimpleAMGTest uses context = " << context->getType() );
    SimpleAMG amg( "AMGTest solver", consoleLogger );
    CriterionPtr minCriterion( new IterationCount( 10 ) );
    amg.setStoppingCriterion( minCriterion );
    amg.setMaxLevels( 2 );
    amg.initialize( coefficients );
    amg.solve( solution, rhs );
    DenseVector<ValueType> diff( solution - refSolution );
    Scalar maxDiff = maxNorm( diff );
    BOOST_CHECK( maxDiff.getValue<ValueType>() < 1e-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( solverTest, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        solverTestMethod<CSRSparseMatrix<ValueType> >( context );
        solverTestMethod<ELLSparseMatrix<ValueType> >( context );
        solverTestMethod<JDSSparseMatrix<ValueType> >( context );
        solverTestMethod<COOSparseMatrix<ValueType> >( context );
        solverTestMethod<DIASparseMatrix<ValueType> >( context );
        solverTestMethod<DenseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE ( testDefaultCriterionSet, ValueType, test_types )
{
    SimpleAMG amg( "SimpleAMG" );
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    const DenseVector<ValueType> rhs( coefficients.getRowDistributionPtr(), 1.0 );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 1.0 );
    DenseVector<ValueType> exactSolution( solution );
    amg.setMaxLevels( 2 );
    amg.initialize( coefficients );
    amg.solve( solution, rhs );
    BOOST_CHECK_EQUAL( amg.getIterationCount(), 1 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    SimpleAMG amg( "SimpleAMG" );
    LAMA_WRITEAT_TEST( amg );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    SimpleAMG amgSolver1( "AMGTestSolver" );
    SolverPtr solverptr = amgSolver1.copy();
    BOOST_CHECK_EQUAL( solverptr->getId(), "AMGTestSolver" );
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

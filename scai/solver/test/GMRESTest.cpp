/**
 * @file GMRESTest.cpp
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
 * @brief Contains the implementation of the class GMRESTest.
 * @author: Malte Förster
 * @date 10.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/test/TestMacros.hpp>

#include <scai/solver/GMRES.hpp>
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

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GMRESTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.GMRESTest" )

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionmethod()
{
    typedef typename MatrixType::MatrixValueType ValueType;
//    LoggerPtr slogger( new CommonLogger(
//        "<GMRES>: ",
//        lama::LogLevel::solverInformation,
//        lama::LoggerWriteBehaviour::toConsoleOnly ) );
    GMRES gmresSolver( "GMRESTestSolver"/*, slogger*/ );
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    DenseVector<ValueType> solution( coefficients.getDistributionPtr(), 1.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getDistributionPtr(), 2.0 );
    DenseVector<ValueType> rhs( coefficients * exactSolution );
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    gmresSolver.setStoppingCriterion( criterion );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    gmresSolver.setPreconditioner( preconditioner );
    gmresSolver.initialize( coefficients );
    gmresSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, gmresSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator(); // default communicator
    testSolveWithPreconditionmethod< CSRSparseMatrix<ValueType> >();
    testSolveWithPreconditionmethod< ELLSparseMatrix<ValueType> >();
    testSolveWithPreconditionmethod< COOSparseMatrix<ValueType> >();
    testSolveWithPreconditionmethod< JDSSparseMatrix<ValueType> >();

// ToDo: these tests do not run with multiple processors

    if ( comm->getSize() <= 1 )
    {
        testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >();
        testSolveWithPreconditionmethod< DenseMatrix<ValueType> >();
    }
}

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionmethod()
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
//    LoggerPtr slogger( new CommonLogger(
//        "<GMRES>: ",
//        lama::LogLevel::solverInformation,
//        lama::LoggerWriteBehaviour::toConsoleOnly ) );
    GMRES gmresSolver( "GMRESTestSolver"/*, slogger */ );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    MatrixType coefficients( helpcoefficients );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );
    const DenseVector<ValueType> rhs( coefficients * exactSolution );
    //initialize
    CriterionPtr criterion( new IterationCount( 10 ) );
    gmresSolver.setStoppingCriterion( criterion );
    gmresSolver.initialize( coefficients );
    gmresSolver.solve( solution, rhs );
    IndexType expectedIterations = 10;
    BOOST_CHECK_EQUAL( expectedIterations, gmresSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator(); // default one
    testSolveWithoutPreconditionmethod< CSRSparseMatrix<ValueType> >();
    testSolveWithoutPreconditionmethod< ELLSparseMatrix<ValueType> >();
    testSolveWithoutPreconditionmethod< JDSSparseMatrix<ValueType> >();
    testSolveWithoutPreconditionmethod< COOSparseMatrix<ValueType> >();

// ToDo: these tests do not run with multiple processors

    if ( comm->getSize() <= 1 )
    {
        testSolveWithoutPreconditionmethod< DIASparseMatrix<ValueType> >();
        testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >();
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testDefaultCriterionSet, ValueType, test_types )
{
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    GMRES gmresSolver( "GMRESTestSolver" );
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> rhs( coefficients.getDistributionPtr(), 2.0 );
    gmresSolver.initialize( coefficients );
    gmresSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( gmresSolver.getIterationCount(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    GMRES gmresSolver( "GMRESTestSolver" );
    LAMA_WRITEAT_TEST( gmresSolver );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

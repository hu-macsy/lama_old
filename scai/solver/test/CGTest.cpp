/**
 * @file CGTest.cpp
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
 * @brief Contains the implementation of the class CGTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 21.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/CG.hpp>
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

#include <scai/hmemo/Context.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CGTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger(
        new CommonLogger( "<CG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    CG cgSolver( "CGTestSolver", slogger );
    BOOST_CHECK_EQUAL( cgSolver.getId(), "CGTestSolver" );
    CG cgSolver2( "CGTestSolver2" );
    BOOST_CHECK_EQUAL( cgSolver2.getId(), "CGTestSolver2" );
    CG cgSolver3( cgSolver2 );
    BOOST_CHECK_EQUAL( cgSolver3.getId(), "CGTestSolver2" );
    BOOST_CHECK( cgSolver3.getPreconditioner() == 0 );
    CG cgSolver4( "cgSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    cgSolver4.setPreconditioner( preconditioner );
    CriterionPtr criterion( new IterationCount( 10 ) );
    cgSolver4.setStoppingCriterion( criterion );
    CG cgSolver5( cgSolver4 );
    BOOST_CHECK_EQUAL( cgSolver5.getId(), cgSolver4.getId() );
    BOOST_CHECK_EQUAL( cgSolver5.getPreconditioner()->getId(), cgSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testSolveWithPreconditionMethod )
{
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<MatrixCreateKeyType> keys;

    Matrix::getCreateValues( keys );

    for ( size_t i = 0; i < keys.size(); ++i )
    {
        if ( keys[i].first == Format::DENSE )
        {
            // does not matter, takes too much time
            continue;
        }

        MatrixPtr coefficients( Matrix::create( keys[i] ) );
        coefficients->setContextPtr( context );

        DistributionPtr dist( new BlockDistribution( coefficients->getNumRows(), comm ) );
        coefficients->redistribute( dist, dist );

        SCAI_LOG_INFO( logger, "CG (with preconditioning) for " << *coefficients << " on " << *context )

        LoggerPtr slogger(
             new CommonLogger( "<CG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly) );
         CG cgSolver( "CGTestSolver", slogger );
         const IndexType N1 = 4;
         const IndexType N2 = 4;
         SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
         CSRSparseMatrix<float> helpcoefficients;
         MatrixCreator<float>::buildPoisson2D( helpcoefficients, 9, N1, N2 );

         // convert to the corresponding matrix type, keep distribution

         *coefficients = helpcoefficients;
         SCAI_LOG_DEBUG( logger, "coefficients matrix = " << *coefficients );

         VectorPtr solution( Vector::getDenseVector( keys[i].second, coefficients->getColDistributionPtr() ) );
         solution->setContextPtr( context );
         *solution = 1.0;

         VectorPtr exactSolution( Vector::getDenseVector( keys[i].second, coefficients->getColDistributionPtr() ) );
         exactSolution->setContextPtr( context );
         *exactSolution = 2.0;

         VectorPtr rhs( Vector::getDenseVector( keys[i].second, coefficients->getRowDistributionPtr() ) ) ;
         rhs->setContextPtr( context );
         *rhs = *coefficients * *exactSolution;

         IndexType expectedIterations = 10;
         CriterionPtr criterion( new IterationCount( expectedIterations ) );
         cgSolver.setStoppingCriterion( criterion );
         SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
         cgSolver.setPreconditioner( preconditioner );
         cgSolver.initialize( *coefficients );
         cgSolver.solve( *solution, *rhs );
         BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );

         VectorPtr diff( Vector::getDenseVector( keys[i].second, coefficients->getColDistributionPtr() ) );
         diff->setContextPtr( context );
         *diff = *solution - *exactSolution;
         Scalar s = l2Norm( *diff );
         SCAI_LOG_INFO( logger,
                        "l2Norm of diff = " << *diff << " = ( solution - exactSolution ) = " << s );
         BOOST_CHECK( s.getValue<float>() < 1E-4 );
    }
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<CSRSparseMatrix<float> > MatrixTypes;

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditionMethod, MatrixType, MatrixTypes )
{
    ContextPtr context = Context::getContextPtr();

    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CG cgSolver( "CGTestSolver" );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "coefficient matrix = " << coefficients );
    coefficients.setContextPtr( context );
    SCAI_LOG_INFO( logger, "CGTest uses context = " << context->getType() );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), 1.0 );
    // Question: should be valid: rhs.getDistribution() == coefficients.getDistribution()
    const DenseVector<ValueType> rhs( coefficients * exactSolution );
    SCAI_LOG_INFO( logger, "rhs = " << rhs );
    //initialize
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "maxNorm of ( solution - exactSolution ) = " << s.getValue<ValueType>() );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testDefaultCriterionSet, ValueType, scai_arithmetic_test_types )
{
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CG cgSolver( "CGTestSolver" );
    CSRSparseMatrix<ValueType> coefficients;
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 2.0 );
    const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 0.0 );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( cgSolver.getIterationCount(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    CG cgSolver( "CGTestSolver" );
    LAMA_WRITEAT_TEST( cgSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    CG cgSolver1( "CGTestSolver" );
    SolverPtr solverptr = cgSolver1.copy();
    BOOST_CHECK_EQUAL( solverptr->getId(), "CGTestSolver" );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

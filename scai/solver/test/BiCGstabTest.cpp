/**
 * @file BiCGstabTest.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCGstabTest.cpp
 * @author lschubert
 * @date 07.08.2013
 * @since 1.1.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/solver/BiCGstab.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/test/TestMacros.hpp>


using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace scai::lama;
using namespace scai::solver;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BiCGstabTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.BiCGstabTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    LoggerPtr slogger( new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );
    BOOST_CHECK_EQUAL( BiCGstabSolver.getId(), "BiCGstabTestSolver" );

    BiCGstab BiCGstabSolver2( "BiCGstabTestSolver2" );
    BOOST_CHECK_EQUAL( BiCGstabSolver2.getId(), "BiCGstabTestSolver2" );

    BiCGstab BiCGstabSolver3( BiCGstabSolver2 );
    BOOST_CHECK_EQUAL( BiCGstabSolver3.getId(), "BiCGstabTestSolver2" );
    BOOST_CHECK( BiCGstabSolver3.getPreconditioner() == 0 );

    BiCGstab BiCGstabSolver4( "BiCGstabSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    BiCGstabSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    BiCGstabSolver4.setStoppingCriterion( criterion );

    BiCGstab BiCGstabSolver5( BiCGstabSolver4 );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getId(), BiCGstabSolver4.getId() );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getPreconditioner()->getId(), BiCGstabSolver4.getPreconditioner()->getId() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testDefaultCriterionSet )
{
    ContextPtr context   = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    typedef double ValueType;
    BiCGstab BiCGstabSolver( "TestBiCGstab" );

    const IndexType N1 = 4;
    const IndexType N2 = 4;

    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );

    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );

    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N1, N2 );

    DenseVector<ValueType> rhs( coefficients.getColDistributionPtr(), 1.0 );
    rhs.setContextPtr( context );

    DenseVector<ValueType> solution( rhs );
    solution.setContextPtr( context );
    solution.redistribute( coefficients.getRowDistributionPtr() );

    BiCGstabSolver.initialize( coefficients );   // Not WORKING

    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK_EQUAL( BiCGstabSolver.getIterationCount(), 1 );
}

/*------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE( testSolveWithIterationCount ) {
    typedef SCAI_TEST_TYPE ValueType;
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    LoggerPtr slogger( new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );

    const IndexType N1 = 40;
    const IndexType N2 = 40;

    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "BiCGstabTest uses context = " << context->getType() );

    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );

    const ValueType solutionInitValue = 1.0;
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), solutionInitValue );
    // TODO: use constructor to set context
    solution.setContextPtr( context );


    DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), solutionInitValue+1.0 );
    // TODO: use constructor to set context
    exactSolution.setContextPtr( context );

    DenseVector<ValueType> rhs( coefficients * exactSolution );

    IndexType maxExpectedIterations = 300;
    CriterionPtr criterion( new IterationCount( maxExpectedIterations ) );
    BiCGstabSolver.setStoppingCriterion( criterion );

    // TODO: this should be tested for ALL solvers
    // Solver has to be initialized before solve is called
    //BOOST_CHECK_THROW ( {BiCGstabSolver.solve( solution, rhs );}, common::Exception );

    BiCGstabSolver.initialize( coefficients );
    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK( maxExpectedIterations == BiCGstabSolver.getIterationCount() );

    DenseVector<ValueType> diff( solution - exactSolution );

    Scalar s                  = maxNorm( diff );
    ValueType realMaxNorm     = s.getValue<ValueType>();
    ValueType expectedMaxNorm = 1E-4;

    SCAI_LOG_INFO( logger, "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << realMaxNorm );

    BOOST_CHECK( realMaxNorm < expectedMaxNorm );


    // Test with preconditioner

    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    BiCGstabSolver.setPreconditioner( preconditioner );

    // TODO: this should be tested for ALL preconditioners
    // Solver has to be initialized AFTER a preconditioner is set
    //BOOST_CHECK_THROW ( {BiCGstabSolver.solve( solution, rhs );}, common::Exception );

    BiCGstabSolver.initialize( coefficients );

    solution = solutionInitValue;

    BiCGstabSolver.solve( solution, rhs );

    BOOST_CHECK( maxExpectedIterations == BiCGstabSolver.getIterationCount() );

    diff = solution - exactSolution;

    s               = maxNorm( diff );
    realMaxNorm     = s.getValue<ValueType>();
    expectedMaxNorm = 1E-4;

    SCAI_LOG_INFO( logger, "maxNorm of diff = " << diff << " = ( solution - exactSolution ) = " << realMaxNorm );

    BOOST_CHECK( realMaxNorm < expectedMaxNorm );

    // TODO: Check that iteration count WITH preconditioner is lower than WITHOUT preconditioner

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    BiCGstab BiCGstabSolver( "BiCGstabTestSolver" );
    LAMA_WRITEAT_TEST( BiCGstabSolver );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    BiCGstab BiCGstabSolver1( "BiCGTestSolver" );

    SolverPtr solverptr = BiCGstabSolver1.copy();

    BOOST_CHECK_EQUAL( solverptr->getId(), "BiCGTestSolver" );
}
/* --------------------------------------------------------------------- */



BOOST_AUTO_TEST_SUITE_END();











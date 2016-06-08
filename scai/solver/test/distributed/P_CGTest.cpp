/**
 * @file P_CGTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains the implementation of the class P_CGTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 27.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_CGTestConfig
{
    P_CGTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
    }

    ~P_CGTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_CGTest, P_CGTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_CGTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithoutPreconditionmethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "testSolveWithoutPreconditionmethod<" << typeid( MatrixType ).name() << " at " << *loc );
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CG cgSolver( "CGTestSolver" );
    SCAI_LOG_INFO( logger, "Solver = " << cgSolver );
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "Poisson2D matrix = " << helpcoefficients );
    MatrixType coefficients( helpcoefficients );
    SCAI_LOG_INFO( logger, "Poisson2D matrix (converted to MatrixType)  = " << helpcoefficients );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContextPtr( loc );
    DenseVector<ValueType> solution( dist, 2.0 );
    const DenseVector<ValueType> exactSolution( dist, 1.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;
    //initialize
    IndexType expectedIterations = 15;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );

    typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

    AbsType sval = s.getValue<AbsType>();

    if ( ! ( sval < 1E-6 ) )
    {
        SCAI_LOG_ERROR( logger, "max norm of diff = " << sval << ", should be < 1E-6 " )
    }

    BOOST_CHECK( sval < 1E-6 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithoutPreconditioning, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    testSolveWithoutPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
    testSolveWithoutPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
    testSolveWithoutPreconditionmethod< DIASparseMatrix<ValueType> >( context );
    testSolveWithoutPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
    testSolveWithoutPreconditionmethod< COOSparseMatrix<ValueType> >( context );
    // @todo: does not work as DenseMatrix = SparseMatrix not supported yet
    // testSolveWithoutPreconditionmethod< DenseMatrix<ValueType> >( context );
}

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveWithPreconditionmethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    SCAI_LOG_INFO( logger, "testSolveWithPreconditionmethod<" << typeid( MatrixType ).name() << "> on " << *loc );
    CG cgSolver( "CGTestSolver" );

    if ( SCAI_LOG_INFO_ON( logger ) )
    {
        LoggerPtr slogger( new CommonLogger(
                               "<CG>: ",
                               scai::solver::LogLevel::solverInformation,
                               scai::solver::LoggerWriteBehaviour::toConsoleOnly ) );
        cgSolver.setLogger( slogger );
    }

    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> csrCoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( csrCoefficients, 9, N1, N2 );
    MatrixType coefficients( csrCoefficients );
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContextPtr( loc );
    DenseVector<ValueType> solution( dist, 1.0 );
    const DenseVector<ValueType> exactSolution( dist, 2.0 );
    DenseVector<ValueType> rhs( dist, 0.0 );
    rhs = coefficients * exactSolution;
    IndexType expectedIterations = 10;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );
    cgSolver.setStoppingCriterion( criterion );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    cgSolver.setPreconditioner( preconditioner );
    SCAI_LOG_INFO( logger, "matrix for CG solver = " << coefficients );
    cgSolver.initialize( coefficients );
    cgSolver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterations, cgSolver.getIterationCount() );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    SCAI_LOG_INFO( logger, "max norm ( solution - exactSolution ) = " << s );

    typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

    if ( s.getValue<AbsType>() >= 1E-6 )
    {
        SCAI_LOG_ERROR( logger, "cgSolver for " << coefficients
                        << ": max norm ( solution - exactSolution ) = " << s );
    }

    BOOST_CHECK( s.getValue<AbsType>() < 1E-6 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolveWithPrecondition, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();
    testSolveWithPreconditionmethod< CSRSparseMatrix<ValueType> >( context );
    testSolveWithPreconditionmethod< ELLSparseMatrix<ValueType> >( context );
    testSolveWithPreconditionmethod< COOSparseMatrix<ValueType> >( context );
    testSolveWithPreconditionmethod< DIASparseMatrix<ValueType> >( context );
    testSolveWithPreconditionmethod< JDSSparseMatrix<ValueType> >( context );
    testSolveWithPreconditionmethod< DenseMatrix<ValueType> >( context );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

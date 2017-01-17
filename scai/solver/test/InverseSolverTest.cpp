/**
 * @file InverseSolverTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Specific tests for the solver class InverseSolver.
 * @author Thomas Brandes, Robin Rehrmann
 * @date 22.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/InverseSolver.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( InverseSolverTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.InverseSolverTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<GMRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    InverseSolver InverseSolverSolver( "InverseSolverSolver", slogger );
    BOOST_CHECK_EQUAL( InverseSolverSolver.getId(), "InverseSolverSolver" );
    InverseSolver InverseSolverSolver2( "InverseSolverSolver2" );
    BOOST_CHECK_EQUAL( InverseSolverSolver2.getId(), "InverseSolverSolver2" );
    InverseSolver InverseSolverSolver3( InverseSolverSolver2 );
    BOOST_CHECK_EQUAL( InverseSolverSolver3.getId(), "InverseSolverSolver2" );
}
// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( InverseTest2, ValueType, scai_numeric_test_types )
{
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get4x4SystemA<ValueType>();
    const IndexType n = 4;
    DenseVector<ValueType> solution( n, 1.0 );
    DenseVector<ValueType> solution2( n, 1.0 );
    std::string s = "DataType";
    InverseSolver inverseSolver( "InverseSolverTest<" + s + "> solver" );
    // DenseMatrix<ValueType> inverse = DenseMatrix<ValueType>( system.coefficients );
    DenseMatrix<ValueType> origin = DenseMatrix<ValueType>( system.coefficients );
    DenseMatrix<ValueType> result = DenseMatrix<ValueType>( system.coefficients );
    inverseSolver.initialize( origin );
    const Matrix& inverse = inverseSolver.getInverse();
    origin.matrixTimesMatrix( result, 1.0, inverse, 0.0, result );

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            Scalar scalar = result.getValue( i, j );

            if ( i == j )
            {
                SCAI_CHECK_CLOSE( 1.0, scalar.getValue<ValueType>(), 1 );
            }
            else
            {
                BOOST_CHECK( scalar.getValue<ValueType>() < 1E-6 );
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------

// copied and adapted from IterativeSolverTest

BOOST_AUTO_TEST_CASE( SolveTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    ContextPtr context = Context::getContextPtr();
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    const IndexType N1 = 10;
    const IndexType N2 = 10;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setContextPtr( context );
    MatrixCreator::buildPoisson2D( coefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "InverseTest uses context = " << context->getType() );
    DistributionPtr rowDist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( coefficients.getNumColumns(), comm ) );
    coefficients.redistribute( rowDist, colDist );
    const ValueType solutionInitValue = 1.0;
    DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), solutionInitValue );
    // TODO: use constructor to set context
    solution.setContextPtr( context );
    DenseVector<ValueType> exactSolution( coefficients.getColDistributionPtr(), solutionInitValue + 1.0 );
    // TODO: use constructor to set context
    exactSolution.setContextPtr( context );
    DenseVector<ValueType> rhs( coefficients * exactSolution );
    IndexType maxExpectedIterations = 3000;
    CriterionPtr criterion( new IterationCount( maxExpectedIterations ) );
    SolverPtr solver ( Solver::create( "InverseSolver", "" ) );
    solver->initialize( coefficients );
    solver->solve( solution, rhs );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s                  = maxNorm( diff );
    ValueType realMaxNorm     = s.getValue<ValueType>();
    ValueType expectedMaxNorm = 1E-4;
    SCAI_LOG_INFO( logger, "maxNorm of diff = " << s << " = ( solution - exactSolution ) = " << realMaxNorm );
    BOOST_CHECK( realMaxNorm < expectedMaxNorm );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

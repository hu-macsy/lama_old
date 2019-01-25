/**
 * @file InverseSolverTest.cpp
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
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

using namespace scai;
using namespace lama;
using namespace solver;
using namespace hmemo;
using namespace dmemo;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( InverseSolverTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.InverseSolverTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<GMRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    InverseSolver<ValueType> solver1( "InverseSolver", slogger );
    BOOST_CHECK_EQUAL( solver1.getId(), "InverseSolver" );
    InverseSolver<ValueType> solver2( "solver2" );
    BOOST_CHECK_EQUAL( solver2.getId(), "solver2" );
    InverseSolver<ValueType> solver3( solver2 );
    BOOST_CHECK_EQUAL( solver3.getId(), "solver2" );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( InverseTest2, ValueType, scai_numeric_test_types )
{
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get4x4SystemA<ValueType>();
    const IndexType n = 4;

    auto solution  = denseVector<ValueType>( n, 1 );
    auto solution2 = denseVector<ValueType>( n, 1 );

    std::string s = "DataType";

    InverseSolver<ValueType> inverseSolver( "InverseSolverTest<" + s + "> solver" );

    // DenseMatrix<ValueType> inverse = DenseMatrix<ValueType>( system.coefficients );

    auto origin = convert<DenseMatrix<ValueType>>( system.coefficients );
    auto result = convert<DenseMatrix<ValueType>>( system.coefficients );

    inverseSolver.initialize( origin );
    const Matrix<ValueType>& inverse = inverseSolver.getInverse();
    origin.matrixTimesMatrix( result, 1.0, inverse, 0.0, result );

    auto eps = common::TypeTraits<ValueType>::small();

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            ValueType scalar = result.getValue( i, j );

            if ( i == j )
            {
                BOOST_CHECK( common::Math::abs( scalar - 1 ) < eps  );
            }
            else
            {
                BOOST_CHECK( common::Math::abs( scalar ) < eps  );
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

    CSRSparseMatrix<ValueType> coefficients( context );

    MatrixCreator::buildPoisson2D( coefficients, 9, N1, N2 );

    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "InverseTest uses context = " << context->getType() );

    auto rowDist = std::make_shared<BlockDistribution>( coefficients.getNumRows(), comm );
    auto colDist = std::make_shared<BlockDistribution>( coefficients.getNumColumns(), comm );

    coefficients.redistribute( rowDist, colDist );

    const ValueType solutionInitValue = 1.0;

    auto solution      = denseVector<ValueType>( colDist, solutionInitValue, context );
    auto exactSolution = denseVector<ValueType>( colDist, solutionInitValue + 1.0, context );
    auto rhs           = denseVectorEval( coefficients * exactSolution, context );

    IndexType maxExpectedIterations = 3000;
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( maxExpectedIterations ) );
    SolverPtr<ValueType> solver ( Solver<ValueType>::getSolver( "InverseSolver" ) );
    solver->initialize( coefficients );
    solver->solve( solution, rhs );
    auto diff = denseVectorEval( solution - exactSolution );

    auto realMaxNorm = maxNorm( diff );                // norm returns RealType<ValueType>
    decltype( realMaxNorm) expectedMaxNorm = 1E-4;

    SCAI_LOG_INFO( logger, "maxNorm of diff ( solution - exactSolution ) = " << realMaxNorm );

    BOOST_CHECK( realMaxNorm < expectedMaxNorm );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

/**
 * @file InverseSolverTest.cpp
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
 * @brief Contains the implementation of the class InverseSolverTest.
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * @since 1.0.0
 **/

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

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

typedef boost::mpl::list<float, double> test_types;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( InverseTest2, ValueType, test_types )
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
                BOOST_CHECK_CLOSE( 1.0, scalar.getValue<ValueType>(), 1 );
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
    MatrixCreator<ValueType>::buildPoisson2D( coefficients, 9, N1, N2 );
    SCAI_LOG_INFO( logger, "coefficients matrix = " << coefficients );
    SCAI_LOG_INFO( logger, "InverseTest uses context = " << context->getType() );

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

    IndexType maxExpectedIterations = 3000;
    CriterionPtr criterion( new IterationCount( maxExpectedIterations ) );

    Solver* solver = Solver::create( "InverseSolver", "" );
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

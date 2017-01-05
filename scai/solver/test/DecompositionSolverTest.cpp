/**
 * @file DecompositionSolverTest.cpp
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Specific test routines for the class DecompositionSolver.
 * @author Lauretta Schubert
 * @date 22.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/DecompositionSolver.hpp>
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

#include <scai/common/TypeTraits.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using namespace scai::common;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( DecompositionSolverTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.DecompositionSolverTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<GMRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    DecompositionSolver DecompositionSolverSolver( "DecompositionSolverSolver", slogger );
    BOOST_CHECK_EQUAL( DecompositionSolverSolver.getId(), "DecompositionSolverSolver" );
    DecompositionSolver DecompositionSolverSolver2( "DecompositionSolverSolver2" );
    BOOST_CHECK_EQUAL( DecompositionSolverSolver2.getId(), "DecompositionSolverSolver2" );
    DecompositionSolver DecompositionSolverSolver3( DecompositionSolverSolver2 );
    BOOST_CHECK_EQUAL( DecompositionSolverSolver3.getId(), "DecompositionSolverSolver2" );
}
// ---------------------------------------------------------------------------------------------------------------

typedef boost::mpl::list<SCAI_NUMERIC_TYPES_EXT_HOST> scai_ext_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( DecompositionTest, ValueType, scai_ext_test_types )
{
    if ( TypeTraits<ValueType>::stype == scalar::LONG_DOUBLE ||
         TypeTraits<ValueType>::stype == scalar::LONG_DOUBLE_COMPLEX )
    {
        // skip because not supported by pardiso or cuSolver
        return;
    }

    if ( TypeTraits<IndexType>::stype != scalar::INT )
    {
        // skip because external solver can only deal with IndexType = int
        return;
    }

    // same test matrix as in CSRUtilsTest
    const IndexType ia[] = { 0, 4, 8, 12, 15 };
    const IndexType ja[] = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3 };
    const ValueType values[] = { 3.0,  4.0, -5.0,  6.0,
                                 6.0,  5.0, -6.0, 5.0,
                                 9.0, -4.0,  2.0, 3.0,
                                       2.0, -3.0, 1.0 };
    const ValueType rhsValues[] = { 39.0, 43.0, 6.0, 13.0 };
    const ValueType solValues[] = { 1.0, 2.0, -2.0, 3.0 };
    const IndexType numRows = 4;
    const IndexType nnz = 15;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() > 1 ) 
    {
        return;   // DecompositionSolver not yet parallel
    }

    DistributionPtr dist( new BlockDistribution( numRows, comm ) );

    ContextPtr context = Context::getContextPtr();

    CSRSparseMatrix<ValueType> matrix;
    matrix.setRawCSRData( dist, dist, nnz, ia, ja, values );
    matrix.setContextPtr( context );

    DenseVector<ValueType> rhs( numRows, rhsValues, context );
    DenseVector<ValueType> solution( numRows, 0.0, context );
    rhs.redistribute( dist );
    solution.redistribute( dist );

    DecompositionSolver solver( "DecompositionSolver" );
    solver.initialize( matrix );
    solver.setContextPtr( context );
    solver.solve( solution, rhs );

    {
        ContextPtr host = Context::getHostPtr();
        ReadAccess<ValueType> rSol( solution.getLocalValues(), host );
        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType x = rSol[i] - solValues[i];
            BOOST_CHECK_SMALL( Math::real( x ), TypeTraits<ValueType>::small() );
            BOOST_CHECK_SMALL( Math::imag( x ), TypeTraits<ValueType>::small() );
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

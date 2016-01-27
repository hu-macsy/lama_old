/**
 * @file P_InverseSolverTest.cpp
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
 * @brief Contains the implementation of the class P_InverseSolverTest.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 27.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/solver/InverseSolver.hpp>
#include <scai/lama/solver/logger/Timer.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( P_InverseSolverTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_InverseSolverTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveMethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    LoggerPtr loggerD(
        new CommonLogger( "<InverseSolver>: ", LogLevel::noLogging,
                          LoggerWriteBehaviour::toConsoleOnly ) );
    InverseSolver inverseSolver( "InverseTestSolver", loggerD );
    SCAI_LOG_DEBUG( logger, "inverseSolver created" )
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    SCAI_LOG_DEBUG( logger, "Poisson2D matrix created" << helpcoefficients )
    MatrixType coefficients( helpcoefficients );
    CommunicatorPtr comm = Communicator::getCommunicator();
    DistributionPtr dist( new BlockDistribution( coefficients.getNumRows(), comm ) );
    coefficients.redistribute( dist, dist );
    coefficients.setContextPtr( loc );
    DenseVector<ValueType> solution( dist, 2.0 );
    const DenseVector<ValueType> exactSolution( dist, 1.0 );
    DenseVector<ValueType> rhs( dist, 1.0 );
    rhs = coefficients * exactSolution;
    SCAI_LOG_INFO( logger, "created all stuff for inverse solver" )
    //initialize
    inverseSolver.initialize( coefficients );
    inverseSolver.setContextPtr( loc );
    inverseSolver.solve( solution, rhs );
    DenseVector<ValueType> diff( solution - exactSolution );
    Scalar s = maxNorm( diff );
    BOOST_CHECK( s.getValue<ValueType>() < 1E-6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveMethod< DIASparseMatrix<ValueType> >( context );
        testSolveMethod< JDSSparseMatrix<ValueType> >( context );
        testSolveMethod< COOSparseMatrix<ValueType> >( context );
        // @todo: DenseMatrix = SparseMatrix not available yet
        // testSolveMethod< DenseMatrix<ValueType> >( context );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

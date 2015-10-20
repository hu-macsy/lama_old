/**
 * @file P_SpecializedJacobi.cpp
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
 * @brief Contains the implementation of the class P_SpecializedJacobi.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 27.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/solver/SpecializedJacobi.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/solver/logger/Timer.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/Communicator.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_SpecializedJacobiTestConfig
{
    P_SpecializedJacobiTestConfig()
    {
        LoggerPtr loggerD(
            new CommonLogger( "<Jacobi>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toConsoleOnly ) );
        mJacobiDouble = new SpecializedJacobi( "SpecializedJacobiTest double solver", loggerD );
        mJacobiFloat = new SpecializedJacobi( "SpecializedJacobiTest float solver", loggerD );
        comm = Communicator::get( "MPI" );
    }

    ~P_SpecializedJacobiTestConfig()
    {
        delete mJacobiDouble;
        delete mJacobiFloat;
        comm = CommunicatorPtr();
    }

    OmegaSolver* mJacobiDouble;
    OmegaSolver* mJacobiFloat;

};

BOOST_FIXTURE_TEST_SUITE( P_SpecializedJacobiTest, P_SpecializedJacobiTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_SpecializedJacobi" );

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void testSolveMethod( ContextPtr loc )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 4;
    const IndexType N2 = 4;
    CSRSparseMatrix<ValueType> helpcoefficients;
    MatrixCreator<ValueType>::buildPoisson2D( helpcoefficients, 9, N1, N2 );
    DistributionPtr dist( new BlockDistribution( helpcoefficients.getNumRows(), comm ) );
    helpcoefficients.redistribute( dist, dist );
    MatrixType coefficients( helpcoefficients );
    coefficients.setContextPtr( loc );
    std::stringstream loggerName;
    loggerName << " <SpecializedJacobi<" << typeid( coefficients ).name() << ">> ";
//    LoggerPtr slogger( new CommonLogger(
//                        loggerName.str(),
//                        LogLevel::solverInformation, //solverInformation, //noLogging,
//                        LoggerWriteBehaviour::toConsoleOnly ) );
    SpecializedJacobi jacobiSolver( "SpecializedJacobiTest"/*, slogger */ );
    DenseVector<ValueType> solution( dist, 1.0 );
    DenseVector<ValueType> exactSolution( solution );
    DenseVector<ValueType> rhs( coefficients * solution );
    SCAI_LOG_INFO( logger, "Matrix for solver: " << coefficients )
    jacobiSolver.initialize( coefficients );
    CriterionPtr criterion( new IterationCount( 40 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    solution = 0.0;
    SCAI_LOG_INFO( logger, "Specialized Jacobi Solver:solve" )
    jacobiSolver.solve( solution, rhs );
    SCAI_LOG_INFO( logger, "l2norm( compute solution - exactSolution )" )
    DenseVector<ValueType> diff( solution - exactSolution );
    L2Norm l2Norm;
    Scalar norm = l2Norm( diff );
    BOOST_CHECK( norm < 1e-1 );
    //bad omega
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        testSolveMethod< CSRSparseMatrix<ValueType> >( context );
        testSolveMethod< ELLSparseMatrix<ValueType> >( context );
        testSolveMethod< JDSSparseMatrix<ValueType> >( context );
        testSolveMethod< COOSparseMatrix<ValueType> >( context );
        testSolveMethod< DIASparseMatrix<ValueType> >( context );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

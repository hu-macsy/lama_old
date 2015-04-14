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

#include <lama/solver/SpecializedJacobi.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/logger/Timer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>

#include <lama/DenseVector.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/Communicator.hpp>

#include <lama/norm/L2Norm.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_SpecializedJacobiTestConfig
{
    P_SpecializedJacobiTestConfig()
    {
        Timer* timerD = new Timer();
        std::auto_ptr<Timer> autoTimerD( timerD );
        Timer* timerF = new Timer();
        std::auto_ptr<Timer> autoTimerF( timerF );
        LoggerPtr loggerD(
            new CommonLogger( "<Jacobi>: ", lama::LogLevel::completeInformation,
                              lama::LoggerWriteBehaviour::toConsoleOnly,
                              std::auto_ptr<Timer>( new Timer() ) ) );
        mJacobiDouble = new SpecializedJacobi( "SpecializedJacobiTest double solver", loggerD );
        mJacobiFloat = new SpecializedJacobi( "SpecializedJacobiTest float solver", loggerD );
        comm = CommunicatorFactory::get( "MPI" );
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

LAMA_LOG_DEF_LOGGER( logger, "Test.P_SpecializedJacobi" );

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
    coefficients.setContext( loc );
    std::stringstream loggerName;
    loggerName << " <SpecializedJacobi<" << typeid( coefficients ).name() << ">> ";
//    LoggerPtr slogger( new CommonLogger(
//                        loggerName.str(),
//                        LogLevel::solverInformation, //solverInformation, //noLogging,
//                        LoggerWriteBehaviour::toConsoleOnly,
//                        std::auto_ptr<Timer>( new Timer() ) ) );
    SpecializedJacobi jacobiSolver( "SpecializedJacobiTest"/*, slogger */ );
    DenseVector<ValueType> solution( dist, 1.0 );
    DenseVector<ValueType> exactSolution( solution );
    DenseVector<ValueType> rhs( coefficients * solution );
    LAMA_LOG_INFO( logger, "Matrix for solver: " << coefficients )
    jacobiSolver.initialize( coefficients );
    CriterionPtr criterion( new IterationCount( 40 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    solution = 0.0;
    LAMA_LOG_INFO( logger, "Specialized Jacobi Solver:solve" )
    jacobiSolver.solve( solution, rhs );
    LAMA_LOG_INFO( logger, "l2norm( compute solution - exactSolution )" )
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

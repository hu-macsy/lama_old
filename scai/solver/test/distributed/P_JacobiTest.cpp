/**
 * @file P_JacobiTest.cpp
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
 * @brief Contains the implementation of the class P_Jacobi.
 * @author Alexander Büchel, Matthias Makulla
 * @date 27.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/Jacobi.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Communicator.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

/* ------------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_JacobiTestConfig
{
    P_JacobiTestConfig()
    {
        LoggerPtr loggerD(
            new CommonLogger( "<Jacobi>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toConsoleOnly ) );
        mJacobiDouble = new Jacobi( "JacobiTest double solver", loggerD );
        mJacobiFloat = new Jacobi( "JacobiTest float solver", loggerD );
        comm = Communicator::getCommunicatorPtr();
    }

    ~P_JacobiTestConfig()
    {
        delete mJacobiDouble;
        delete mJacobiFloat;
        comm = CommunicatorPtr();
    }

    OmegaSolver* mJacobiDouble;
    OmegaSolver* mJacobiFloat;

};

BOOST_FIXTURE_TEST_SUITE( P_JacobiTest, P_JacobiTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_Jacobi" );

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
    loggerName << " <Jacobi<" << typeid( coefficients ).name() << ">> ";
//    LoggerPtr slogger( new CommonLogger(
//                        loggerName.str(),
//                        LogLevel::solverInformation, //solverInformation, //noLogging,
//                        LoggerWriteBehaviour::toConsoleOnly ) );
    Jacobi jacobiSolver( "JacobiTest"/*, slogger */ );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( testSolve, ValueType, scai_arithmetic_test_types )
{
    if ( scai::common::isComplex( scai::common::TypeTraits<ValueType>::stype ) )
    {
        return;    // do not test for complex types
    }

    ContextPtr context = Context::getContextPtr();

    testSolveMethod< CSRSparseMatrix<ValueType> >( context );
    testSolveMethod< ELLSparseMatrix<ValueType> >( context );
    testSolveMethod< JDSSparseMatrix<ValueType> >( context );
    testSolveMethod< COOSparseMatrix<ValueType> >( context );
    testSolveMethod< DIASparseMatrix<ValueType> >( context );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

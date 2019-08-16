/**
 * @file SimpleAMGTest.cpp
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
 * @brief Specific test routines for the solver class SimpleAMG.
 * @author Malte Foerster
 * @date 22.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/FileLogger.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai;

using namespace lama;
using namespace solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( SimpleAMGTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SimpleAMGTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<SimpleAMG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    SimpleAMG<ValueType> cgSolver( "SimpleAMGTestSolver", slogger );
    BOOST_CHECK_EQUAL( cgSolver.getId(), "SimpleAMGTestSolver" );
    SimpleAMG<ValueType> cgSolver2( "SimpleAMGTestSolver2" );
    BOOST_CHECK_EQUAL( cgSolver2.getId(), "SimpleAMGTestSolver2" );
    SimpleAMG<ValueType> cgSolver3( cgSolver2 );
    BOOST_CHECK_EQUAL( cgSolver3.getId(), "SimpleAMGTestSolver2" );
    BOOST_CHECK( cgSolver3.getPreconditioner() == 0 );
    SimpleAMG<ValueType> cgSolver4( "cgSolver4" );
    SolverPtr<ValueType> preconditioner( new TrivialPreconditioner<ValueType>( "Trivial preconditioner" ) );
    cgSolver4.setPreconditioner( preconditioner );
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 10 ) );
    cgSolver4.setStoppingCriterion( criterion );
    SimpleAMG<ValueType> cgSolver5( cgSolver4 );
    BOOST_CHECK_EQUAL( cgSolver5.getId(), cgSolver4.getId() ); 
    BOOST_CHECK_EQUAL( cgSolver5.getPreconditioner()->getId(), cgSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE ( SetterTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    const IndexType N = 40;
    scai::lama::CSRSparseMatrix<ValueType> coefficients;
    scai::lama::MatrixCreator::buildPoisson2D( coefficients, 5, N, N );
    SimpleAMG<ValueType> SimpleAMGSolver( "SimpleAMGSolver" );

    const IndexType hostOnlyLevel = 13;
    SimpleAMGSolver.setHostOnlyLevel( hostOnlyLevel );

    const IndexType hostOnlyVars = 1150;
    SimpleAMGSolver.setHostOnlyVars( hostOnlyVars );

    const IndexType maxLevels = 25;
    SimpleAMGSolver.setMaxLevels( maxLevels );

    const IndexType minVarsCoarseLevel = 100;
    SimpleAMGSolver.setMinVarsCoarseLevel( minVarsCoarseLevel );

    const IndexType replicatedLevel = 8;
    SimpleAMGSolver.setReplicatedLevel( replicatedLevel );

    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSetup().getHostOnlyLevel(), hostOnlyLevel );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSetup().getHostOnlyVars(), hostOnlyVars );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSetup().getMaxLevels(), maxLevels );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSetup().getMinVarsCoarseLevel(), minVarsCoarseLevel );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSetup().getReplicatedLevel(), replicatedLevel );

    SolverPtr<ValueType> cgSolver ( new CG<ValueType> ( "CGCoarseLevelSolver" ) );
    SimpleAMGSolver.setCoarseLevelSolver( cgSolver );

    SolverPtr<ValueType> jacobiSolver ( new Jacobi<ValueType> ( "JacobiSmoother" ) );
    SimpleAMGSolver.setSmoother( jacobiSolver );

    SimpleAMGSolver.initialize( coefficients ); 

    BOOST_CHECK_EQUAL( SimpleAMGSolver.getCoarseLevelSolver().getId(), cgSolver->getId() );

    for ( IndexType i = 0; i < SimpleAMGSolver.getNumLevels() - 1; ++i )
    {
        BOOST_CHECK_EQUAL( SimpleAMGSolver.getSmoother( i ).getId(), jacobiSolver->getId() );
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE ( SolveTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 40;

    scai::lama::CSRSparseMatrix<ValueType> coefficients;
    scai::lama::MatrixCreator::buildPoisson2D( coefficients, 5, N, N );

    const dmemo::Communicator& comm = coefficients.getRowDistribution().getCommunicator();

    SCAI_LOG_INFO( logger, "Solve matrix: " << coefficients )

    std::string loggerFileName = "logger%r.solve";

    // use class PartitionIO to get individual file names for each processor

    std::string pLoggerFileName = loggerFileName;
    bool isPartitioned;
    PartitionIO::getPartitionFileName( pLoggerFileName, isPartitioned, comm );

    std::shared_ptr<Timer> timer( new Timer() );

    LoggerPtr slogger( new CommonLogger( "<SimpleAMG>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toFileOnly, pLoggerFileName, timer ) );
    SimpleAMG<ValueType> simpleAMGSolver( "SimpleAMGSolver", slogger );

    simpleAMGSolver.initialize( coefficients );

    auto rhs = denseVector<ValueType>( coefficients.getRowDistributionPtr(), 1 );
    auto x = denseVector<ValueType>( coefficients.getColDistributionPtr(), 0 );

    simpleAMGSolver.solve( x, rhs );

    // check that timing works

    if ( simpleAMGSolver.getNumLevels() > 1 )
    {
        BOOST_CHECK( simpleAMGSolver.getAverageSmootherTime() > 0.0 );
        BOOST_CHECK( simpleAMGSolver.getAverageTransferTime() > 0.0 );
        BOOST_CHECK( simpleAMGSolver.getAverageResidualTime() > 0.0 );
    }

    // Now delete the logger file, close it before

    FileLogger::getFileLogger().closeLogFile();

    int rc = PartitionIO::removeFile( loggerFileName, comm );
    BOOST_CHECK_EQUAL( 0, rc );
    BOOST_CHECK( !PartitionIO::fileExists( loggerFileName, comm ) );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

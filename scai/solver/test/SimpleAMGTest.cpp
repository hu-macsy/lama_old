/**
 * @file SimpleAMGTest.cpp
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

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<SimpleAMG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    SimpleAMG SimpleAMGSolver( "SimpleAMGSolver", slogger );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getId(), "SimpleAMGSolver" );
    SimpleAMG SimpleAMGSolver2( "SimpleAMGSolver2" );
    BOOST_CHECK_EQUAL( SimpleAMGSolver2.getId(), "SimpleAMGSolver2" );
    SimpleAMG SimpleAMGSolver3( SimpleAMGSolver2 );
    BOOST_CHECK_EQUAL( SimpleAMGSolver3.getId(), "SimpleAMGSolver2" );
    BOOST_CHECK( SimpleAMGSolver3.getPreconditioner() == 0 );
    SimpleAMG SimpleAMGSolver4( "SimpleAMGSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    SimpleAMGSolver4.setPreconditioner( preconditioner );
    CriterionPtr criterion( new IterationCount( 10 ) );
    SimpleAMGSolver4.setStoppingCriterion( criterion );
    SimpleAMG SimpleAMGSolver5( SimpleAMGSolver4 );
    BOOST_CHECK_EQUAL( SimpleAMGSolver5.getId(), SimpleAMGSolver4.getId() );
    BOOST_CHECK_EQUAL( SimpleAMGSolver5.getPreconditioner()->getId(), SimpleAMGSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

//TODO:
BOOST_AUTO_TEST_CASE ( SetterTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    const IndexType N = 40;
    scai::lama::CSRSparseMatrix<ValueType> coefficients;
    scai::lama::MatrixCreator::buildPoisson2D( coefficients, 5, N, N );
    SimpleAMG SimpleAMGSolver( "SimpleAMGSolver" );
    /* cant not be tested: not getter
    SimpleAMGSolver.setHostOnlyLevel( IndexType hostOnlyLevel );
    SimpleAMGSolver.setHostOnlyVars( IndexType hostOnlyVars );
    SimpleAMGSolver.setReplicatedLevel( IndexType replicatedLevel );
    SimpleAMGSolver.setMaxLevels( unsigned int levels );
    SimpleAMGSolver.setMinVarsCoarseLevel( unsigned int vars );*/
    SolverPtr cgSolver ( new CG ( "CGCoarseLevelSolver" ) );
    scai::hmemo::ContextPtr context = scai::hmemo::Context::getContextPtr();
    // does not work because coarselevelsolver in SingleGridSetup gets overridden by smoother
    //SimpleAMGSolver.setCoarseLevelSolver( cgSolver );
    //SimpleAMGSolver.initialize( coefficients ); // solver needs to be initialized to have a AMGSetup to pass coarselevelsolver to setup
    //BOOST_CHECK_EQUAL( SimpleAMGSolver.getCoarseLevelSolver().getId(), cgSolver->getId() );
    SimpleAMGSolver.setSmoother( cgSolver );
    SimpleAMGSolver.initialize( coefficients ); // solver needs to be initialized to have a AMGSetup to pass coarselevelsolver to setup
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSmoother( 0 ).getId(), cgSolver->getId() );
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

    common::shared_ptr<Timer> timer( new Timer() );

    LoggerPtr slogger( new CommonLogger( "<SimpleAMG>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toFileOnly, pLoggerFileName, timer ) );
    SimpleAMG simpleAMGSolver( "SimpleAMGSolver", slogger );

    simpleAMGSolver.initialize( coefficients );

    for ( unsigned int i = 0; i < simpleAMGSolver.getNumLevels(); ++i )
    {
        const _Matrix& galerkin    = simpleAMGSolver.getGalerkin( i );
        const _Matrix& restriction = simpleAMGSolver.getRestriction( i );
        const _Matrix& interpol    = simpleAMGSolver.getInterpolation( i );
        const Solver& smoother    = simpleAMGSolver.getSmoother( i );

        if ( i == 0 )
        {
            BOOST_CHECK_EQUAL( galerkin.getNumRows(), coefficients.getNumRows() );
        }

        BOOST_CHECK_EQUAL( restriction.getNumRows(), interpol.getNumColumns() );
        BOOST_CHECK_EQUAL( interpol.getNumRows(), restriction.getNumColumns() );

        // smoother must be an iterative solver

        BOOST_CHECK( dynamic_cast<const IterativeSolver*>( &smoother ) );
    }

    DenseVector<ValueType> rhs( coefficients.getRowDistributionPtr() );
    DenseVector<ValueType> x( coefficients.getColDistributionPtr() );

    rhs = ValueType( 1 );
    x = ValueType( 0 );

    simpleAMGSolver.solve( x, rhs );

    // check that timing works

    BOOST_CHECK( simpleAMGSolver.getAverageSmootherTime() > 0.0 );
    BOOST_CHECK( simpleAMGSolver.getAverageTransferTime() > 0.0 );
    BOOST_CHECK( simpleAMGSolver.getAverageResidualTime() > 0.0 );

    // Now delete the logger file, close it before

    FileLogger::getFileLogger().closeLogFile();

    int rc = PartitionIO::removeFile( loggerFileName, comm );
    BOOST_CHECK_EQUAL( 0, rc );
    BOOST_CHECK( !PartitionIO::fileExists( loggerFileName, comm ) );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

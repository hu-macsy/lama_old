/**
 * @file SimpleAMGTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Contains the implementation of the class SimpleAMGTest.cpp
 * @author Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/CG.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

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
    scai::lama::MatrixCreator<ValueType>::buildPoisson2D( coefficients, 5, N, N );

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
    SimpleAMGSolver.setSmootherContext( context );
    SimpleAMGSolver.initialize( coefficients ); // solver needs to be initialized to have a AMGSetup to pass coarselevelsolver to setup
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSmoother( 0 ).getId(), cgSolver->getId() );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getSmoother( 0 ).getContextPtr()->getType(), context->getType() );
}

// ---------------------------------------------------------------------------------------------------------------

// GetterTest makes no sense, because it test the AMGSetup, which can be SingleGridSetup or other, with different defaults
/*    const lama::Matrix& getGalerkin( unsigned int level );
    const lama::Matrix& getRestriction( unsigned int level );
    const lama::Matrix& getInterpolation( unsigned int level );

    lama::Vector& getSolutionVector( unsigned int level );
    lama::Vector& getRhsVector( unsigned int level );

    Solver& getSmoother( unsigned int level );
    Solver& getCoarseLevelSolver();

    double getAverageSmootherTime() const;
    double getAverageTransferTime() const;
    double getAverageResidualTime() const;
*/

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

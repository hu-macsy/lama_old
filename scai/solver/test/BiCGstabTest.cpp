/**
 * @file BiCGstabTest.cpp
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
 * @brief BiCGstabTest.cpp
 * @author lschubert
 * @date 07.08.2013
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/BiCGstab.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/test/TestMacros.hpp>


using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( BiCGstabTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.BiCGstabTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<BiCGstab>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    BiCGstab BiCGstabSolver( "BiCGstabTestSolver", slogger );
    BOOST_CHECK_EQUAL( BiCGstabSolver.getId(), "BiCGstabTestSolver" );

    BiCGstab BiCGstabSolver2( "BiCGstabTestSolver2" );
    BOOST_CHECK_EQUAL( BiCGstabSolver2.getId(), "BiCGstabTestSolver2" );

    BiCGstab BiCGstabSolver3( BiCGstabSolver2 );
    BOOST_CHECK_EQUAL( BiCGstabSolver3.getId(), "BiCGstabTestSolver2" );
    BOOST_CHECK( BiCGstabSolver3.getPreconditioner() == 0 );

    BiCGstab BiCGstabSolver4( "BiCGstabSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    BiCGstabSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    BiCGstabSolver4.setStoppingCriterion( criterion );

    BiCGstab BiCGstabSolver5( BiCGstabSolver4 );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getId(), BiCGstabSolver4.getId() );
    BOOST_CHECK_EQUAL( BiCGstabSolver5.getPreconditioner()->getId(), BiCGstabSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

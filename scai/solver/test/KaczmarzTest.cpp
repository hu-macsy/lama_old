/**
 * @file KaczmarzTest.cpp
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
 * @brief Contains the implementation of the class KaczmarzTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 21.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/Kaczmarz.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( KaczmarzTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.KaczmarzTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<Kaczmarz>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    Kaczmarz cgSolver( "KaczmarzTestSolver", slogger );
    BOOST_CHECK_EQUAL( cgSolver.getId(), "KaczmarzTestSolver" );
    Kaczmarz cgSolver2( "KaczmarzTestSolver2" );
    BOOST_CHECK_EQUAL( cgSolver2.getId(), "KaczmarzTestSolver2" );
    Kaczmarz cgSolver3( cgSolver2 );
    BOOST_CHECK_EQUAL( cgSolver3.getId(), "KaczmarzTestSolver2" );
    BOOST_CHECK( cgSolver3.getPreconditioner() == 0 );
    Kaczmarz cgSolver4( "cgSolver4" );
    CriterionPtr criterion( new IterationCount( 10 ) );
    cgSolver4.setStoppingCriterion( criterion );
    Kaczmarz cgSolver5( cgSolver4 );
    BOOST_CHECK_EQUAL( cgSolver5.getId(), cgSolver4.getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
/**
 * @file KaczmarzTest.cpp
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
 * @brief Specific tests for the solver class Kaczmarz.
 * @author Thomas Brandes
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

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<Kaczmarz>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    Kaczmarz<ValueType> solver( "KaczmarzTestSolver", slogger );
    BOOST_CHECK_EQUAL( solver.getId(), "KaczmarzTestSolver" );
    Kaczmarz<ValueType> solver2( "KaczmarzTestSolver2" );
    BOOST_CHECK_EQUAL( solver2.getId(), "KaczmarzTestSolver2" );
    Kaczmarz<ValueType> solver3( solver2 );
    BOOST_CHECK_EQUAL( solver3.getId(), "KaczmarzTestSolver2" );
    BOOST_CHECK( solver3.getPreconditioner() == 0 );
    Kaczmarz<ValueType> solver4( "solver4" );
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 10 ) );
    solver4.setStoppingCriterion( criterion );
    Kaczmarz<ValueType> solver5( solver4 );
    BOOST_CHECK_EQUAL( solver5.getId(), solver4.getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

/**
 * @file CGSTest.cpp
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
 * @brief Specific tests for the solver class CGS.
 * @author David Schissler
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/CGS.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( CGSTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGSTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<CGS>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    CGS<ValueType> solver( "CGSTestSolver", slogger );
    BOOST_CHECK_EQUAL( solver.getId(), "CGSTestSolver" );
    CGS<ValueType> solver2( "CGSTestSolver2" );
    BOOST_CHECK_EQUAL( solver2.getId(), "CGSTestSolver2" );
    CGS<ValueType> solver3( solver2 );
    BOOST_CHECK_EQUAL( solver3.getId(), "CGSTestSolver2" );
    BOOST_CHECK( solver3.getPreconditioner() == 0 );
    CGS<ValueType> solver4( "solver4" );
    SolverPtr<ValueType> preconditioner( new TrivialPreconditioner<ValueType>( "Trivial preconditioner" ) );
    solver4.setPreconditioner( preconditioner );
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 10 ) );
    solver4.setStoppingCriterion( criterion );
    CGS<ValueType> solver5( solver4 );
    BOOST_CHECK_EQUAL( solver5.getId(), solver4.getId() );
    BOOST_CHECK_EQUAL( solver5.getPreconditioner()->getId(), solver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

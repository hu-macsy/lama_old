/**
 * @file BiCGTest.cpp
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
 * @brief Speficic tests for the solver class BiCG.
 * @author Lauretta Schubert
 * @date 04.07.2013
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/BiCG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( BiCGTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.BiCGTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, ValueType, scai_numeric_test_types )
{
    LoggerPtr slogger( new CommonLogger( "<BiCG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    BiCG<ValueType> solver( "BiCGTestSolver", slogger );
    BOOST_CHECK_EQUAL( solver.getId(), "BiCGTestSolver" );
    BiCG<ValueType> solver2( "BiCGTestSolver2" );
    BOOST_CHECK_EQUAL( solver2.getId(), "BiCGTestSolver2" );
    BiCG<ValueType> solver3( solver2 );
    BOOST_CHECK_EQUAL( solver3.getId(), "BiCGTestSolver2" );
    BOOST_CHECK( solver3.getPreconditioner() == 0 );
    BiCG<ValueType> solver4( "solver4" );
    SolverPtr<ValueType> preconditioner( new TrivialPreconditioner<ValueType>( "Trivial preconditioner" ) );
    solver4.setPreconditioner( preconditioner );
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 10 ) );
    solver4.setStoppingCriterion( criterion );
    BiCG<ValueType> solver5( solver4 );
    BOOST_CHECK_EQUAL( solver5.getId(), solver4.getId() );
    BOOST_CHECK_EQUAL( solver5.getPreconditioner()->getId(), solver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

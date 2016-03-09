/**
 * @file BiCGstabTest.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCGstabTest.cpp
 * @author lschubert
 * @date 07.08.2013
 * @since 1.1.0
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

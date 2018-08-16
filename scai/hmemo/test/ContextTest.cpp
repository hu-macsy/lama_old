/**
 * @file hmemo/test/ContextTest.cpp
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
 * @brief ToDo: Missing description in ./hmemo/test/ContextTest.cpp
 * @author Thomas Brandes
 * @date 08.07.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/ContextAccess.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include "MockContext.hpp"

using namespace scai;
using namespace hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ContextTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ContextTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( hostContextText )
{
    // make sure that host context is always available
    BOOST_CHECK( Context::canCreate( common::ContextType::Host ) );
    ContextPtr host = Context::create( common::ContextType::Host, -1 );
    BOOST_CHECK( host.get() );
    BOOST_CHECK( host == Context::getHostPtr() );
    SCAI_LOG_INFO( logger, "host context = " << *host )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextText )
{
    // some basic test to get the specified context
    ContextPtr ctx = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "context = " << *ctx << ", type = " << ctx->getType() )
    BOOST_CHECK( ctx.get() );
    BOOST_CHECK( Context::hasContext( ctx->getType() ) );
    // context of this type should be available directly
    ContextPtr ctx1 = Context::getContextPtr( ctx->getType() );
    // equality of devices checked by pointer equality
    BOOST_CHECK_EQUAL( ctx.get(), ctx1.get() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( currentContextTest )
{
    ContextPtr ctx = Context::getContextPtr();

    BOOST_CHECK( Context::getCurrentContext() == NULL );

    {
        SCAI_CONTEXT_ACCESS( ctx )
        BOOST_CHECK( Context::getCurrentContext() == ctx.get() );
    }

    BOOST_CHECK( Context::getCurrentContext() == NULL );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( useContextTest )
{
    ContextPtr userContext  = Context::getContextPtr( common::ContextType::UserContext, 1 );
    ContextPtr hostContext  = Context::getContextPtr( common::ContextType::Host );
    ContextPtr testContext  = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "testContext = " << *testContext << ", userContext = " << *userContext );
    const IndexType N = 7;
    HArray<double> X( N, 5.0 );
    // take ownership on userContext
    {
        WriteAccess<double> write( X, userContext );
    }
    // take ownership on testContext
    {
        WriteAccess<double> write( X, testContext );
    }
    // read @ userContext: valid data is transfered from testContext to here
    ReadAccess<double> read( X, userContext );
    const double* vals = read.get();

    for ( int i = 0; i < 7; ++i )
    {
        SCAI_ASSERT_EQUAL( vals[i], 5.0, "check" )
    }

    HArray<double> Y( X );
    Y.clear();
    Y.purge();
    HArray<float> v ( 4, 1.0f );
    {
        ReadAccess<float> read( v, userContext );
        WriteAccess<float> write( v, userContext );
    }

    // read and write access at same time not possible

    if ( userContext.get() != testContext.get() )
    {
        SCAI_LOG_DEBUG( logger, "exception must be thrown for read/write on different devices" )
        BOOST_CHECK_THROW(
        {
            ReadAccess<float> read( v, userContext );
            WriteAccess<float> write( v, testContext );
        },
        common::Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();


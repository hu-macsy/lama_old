/**
 * @file hmemo/test/ContextStackTest.cpp
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
 * @brief Some test for using scope of context via a context stack.
 * @author Thomas Brandes
 * @date 30.04.2019
 */
#include <boost/test/unit_test.hpp>

#include <scai/hmemo/ContextStack.hpp>
#include <scai/common/SCAITypes.hpp>
#include "MockContext.hpp"

using namespace scai;
using namespace hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ContextStackTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ContextStackTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scopeText )
{
    ContextPtr mainContext  = Context::getContextPtr();
    ContextPtr userContext  = Context::getContextPtr( common::ContextType::UserContext, 1 );
    ContextPtr hostContext  = Context::getContextPtr( common::ContextType::Host );

    BOOST_CHECK( *userContext != *hostContext );

    SCAI_LOG_INFO( logger, "actual context: " << *Context::getContextPtr() )

    {
        SCAI_HMEMO_CONTEXT( userContext );
    
        SCAI_LOG_INFO( logger, "actual context: " << *Context::getContextPtr() )

        BOOST_CHECK_EQUAL( *userContext, *Context::getContextPtr() );
 
        {
            SCAI_HMEMO_CONTEXT( hostContext );
        
            SCAI_LOG_INFO( logger, "actual context: " << *Context::getContextPtr() )

            BOOST_CHECK_EQUAL( *hostContext, *Context::getContextPtr() );
        }

        SCAI_LOG_INFO( logger, "actual context: " << *Context::getContextPtr() )

        BOOST_CHECK_EQUAL( *userContext, *Context::getContextPtr() );
    }

    SCAI_LOG_INFO( logger, "actual context: " << *Context::getContextPtr() )

    BOOST_CHECK_EQUAL( *mainContext, *Context::getContextPtr() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();


/**
 * @file hmemo/test/MemoryTest.cpp
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
 * @brief ToDo: Missing description in ./hmemo/test/MemoryTest.cpp
 * @author Thomas Brandes
 * @date 24.02.2016
 */
#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/exception/MemoryException.hpp>

BOOST_AUTO_TEST_SUITE( MemoryTest )

using namespace scai;
using namespace hmemo;

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MemoryTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DeviceMemoryTest )
{
    ContextPtr ctx = Context::getContextPtr();  // actual test context
    MemoryPtr  mem = ctx->getMemoryPtr();
    // make sure that it is not a NULL memory
    BOOST_REQUIRE( mem.get() != NULL );
    SCAI_LOG_INFO( logger, "Test of device memory " << *mem )
    // important: memory of context can be used with context
    BOOST_CHECK( ctx->canUseMemory( *mem ) );
    size_t N = 100;
    void* data = mem->allocate( N );
    BOOST_ASSERT( data != NULL );
    mem->free( data, N );

    size_t NM = mem->maxAllocatedBytes() + 1;
    data = mem->allocate( NM );
    mem->free( data, NM );

    BOOST_ASSERT( mem->maxAllocatedBytes() >= NM );

    // Make sure that too much memory allocation throws an exception
    size_t MAX_N = std::numeric_limits<size_t>::max();
    BOOST_CHECK_THROW(
    {
        data = mem->allocate( MAX_N );
    }
    , MemoryException );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( HostMemoryTest )
{
    ContextPtr ctx = Context::getContextPtr();  // actual test context
    MemoryPtr  mem = ctx->getHostMemoryPtr();
    // make sure that it is not a NULL memory
    BOOST_REQUIRE( mem.get() != NULL );
    SCAI_LOG_INFO( logger, "Test of host memory " << *mem << " for " << *ctx )
    // important: host memory of context can be used on host
    ContextPtr host = Context::getHostPtr();  // host context
    BOOST_CHECK( host->canUseMemory( *mem ) );
    size_t N = 100;
    void* data = mem->allocate( N );
    BOOST_ASSERT( data != NULL );
    mem->free( data, N );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyTest )
{
    ContextPtr ctx  = Context::getContextPtr();  // actual test context
    MemoryPtr  cmem = ctx->getMemoryPtr();
    MemoryPtr  hmem = ctx->getHostMemoryPtr();   // can be usual host or pinned memory
    SCAI_LOG_INFO( logger, "CopyTest memory " << *cmem << " to/from " << *hmem )
    // each memory should be able to copy to and copy from Host
    BOOST_REQUIRE( cmem->canCopyTo( *hmem ) );
    BOOST_REQUIRE( cmem->canCopyFrom( *hmem ) );
    const size_t N = 10;
    char* cdata = reinterpret_cast<char*>( cmem->allocate( N ) );
    char* hdata = reinterpret_cast<char*>( hmem->allocate( N ) );
    // Init host memory with 0
    const char VAL = 13;  // should not be 0
    memset( hdata, VAL, N );
    cmem->memcpyFrom( cdata, *hmem, hdata, N );
    memset( hdata, 0, N );
    cmem->memcpyTo( *hmem, hdata, cdata, N );

    for ( size_t i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( hdata[i], VAL );
    }

    hmem->free( hdata, N );
    cmem->free( cdata, N );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

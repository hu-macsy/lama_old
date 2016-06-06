/**
 * @file MICContextTest.cpp
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
 * @brief Extensive testing of MICContext.
 * @author Thomas Brandes
 * @date 13.09.2015
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/hmemo.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

using namespace scai;

using namespace common;
using namespace hmemo;

typedef boost::mpl::list<double, float> test_types;

BOOST_AUTO_TEST_SUITE( MICContextTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.MICContextTest" );

template<typename ValueType>
static void scal( int n, ValueType alpha, ValueType* array )
{
    void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n, alpha )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );

        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] *= alpha;
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    // Test will take the default MIC device
    ContextPtr micContext1 = Context::getContextPtr( context::MIC );
    ContextPtr micContext2 = Context::getContextPtr( context::MIC );
    // Two queries for the same context should deliver same pointer
    BOOST_CHECK( micContext1.get() == micContext2.get() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getMemoryTest )
{
    // Test will take the default MIC device
    ContextPtr micContext = Context::getContextPtr( context::MIC );
    MemoryPtr micMemory = micContext->getMemoryPtr();
    BOOST_CHECK( micMemory.get() != NULL );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    HArray<int> ctxArray; // default, not allocated at all
    ContextPtr micContext = Context::getContextPtr( context::MIC );
    {
        WriteAccess<int> array( ctxArray, micContext );
        array.resize( 10 );
        // Destructors will be called for array, ctxArray, micContext
    }
    BOOST_CHECK_EQUAL( 10, ctxArray.size() );
    BOOST_CHECK( ctxArray.capacity( micContext ) >= 10 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( useTest )
{
    ContextPtr micContext = Context::getContextPtr( context::MIC );

    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    HArray<float> vector( n, value );
    {
        WriteAccess<float> micV( vector, micContext );
        SCAI_CONTEXT_ACCESS( micContext );
        scal( n, alpha, micV.get() );
    }
    ReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV.get()[i] );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();


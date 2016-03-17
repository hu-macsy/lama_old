/**
 * @file AllStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Test cases applied to each storage class, i.e. test (virtual) methods of _MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;
using utilskernel::LArray;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( AllStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.AllStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    Storages allMatrixStorages;    // is created by factory

    size_t nFormats = Format::UNDEFINED;
    size_t nTypes   = ARITHMETIC_HOST_TYPE_CNT;

    SCAI_LOG_INFO( logger, "Test all storages of factory to be empty, #storages = " << allMatrixStorages.size() )

    BOOST_CHECK_EQUAL( nTypes * nFormats, allMatrixStorages.size() );

    for ( size_t i = 0; i < allMatrixStorages.size(); ++i )
    {
        _MatrixStorage& storage = *allMatrixStorages[i];

        BOOST_CHECK_EQUAL( 0, storage.getNumRows() );
        BOOST_CHECK_EQUAL( 0, storage.getNumColumns() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setIdentityTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    const IndexType N = 15;        // size of the matrix storage

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for setIdentity" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];

        storage.setIdentity(  N );

        SCAI_LOG_DEBUG( logger, "Identity matrix, N = " << N << " : " << storage )

        BOOST_REQUIRE_EQUAL( N, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumColumns() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumValues() );
        
        LArray<double> row;

        for ( IndexType i = 0; i < N; ++i )
        {
            storage.getRow( row, i );
            hmemo::ReadAccess<double> rRow( row );
        
            for ( IndexType j = 0; j < N; ++j )
            {
                if ( i == j )
                {
                    BOOST_CHECK_EQUAL( 1.0, rRow[j] );
                }
                else
                {
                    BOOST_CHECK_EQUAL( 0.0, rRow[j] );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();

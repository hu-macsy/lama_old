/**
 * @file TypedStorageTest.cpp
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
 * @brief Test cases applied to each typed storage class, i.e. test (virtual) methods of MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] = { 6, 0, 0, 4, 7, 0, 0, 0, 0, 0, -9.3f, 4, 2, 5, 0, 3 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( TypedStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.TypedStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( factoryTest, ValueType, scai_arithmetic_test_types )
{
    TypedStorages<ValueType> allMatrixStorages;    // is created by factory

    size_t nFormats = Format::UNDEFINED;

    SCAI_LOG_INFO( logger, "factoryTest<" << common::TypeTraits<ValueType>::id() << "> : " 
                            << allMatrixStorages.size() << " storages"                        )

    BOOST_CHECK_EQUAL( nFormats, allMatrixStorages.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( normTest, ValueType, scai_arithmetic_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    SCAI_LOG_INFO( logger, "normTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_DEBUG( logger, "normTest, storage = " << storage << " @ " << *storage.getContextPtr() )

        setDenseData( storage );

        ValueType maxNorm = storage.maxNorm();

        ValueType expected = 9.3f; // maximal absolute value

        SCAI_CHECK_CLOSE( maxNorm, expected, 1 );
    }        
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_arithmetic_test_types )
{
    DenseStorage<ValueType> tmp;

    setDenseData( tmp );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_DEBUG( logger, "scaleTest, storage = " << storage << " @ " << *storage.getContextPtr() )

        setDenseData( storage );

        storage.scale( 2.0 );  // should be executed on context

        for ( IndexType i = 0; i < storage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < storage.getNumColumns(); ++j )
            {
                BOOST_CHECK_EQUAL( 2.0 * tmp.getValue( i, j ), storage.getValue( i, j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */


BOOST_AUTO_TEST_SUITE_END();

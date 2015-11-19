/**
 * @file MatrixStorageConversionTest.cpp
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
 * @brief Test conversions between different matrix storage formats
 * @author Thomas Brandes
 * @date 15.07.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/SparseAssemblyStorage.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

template<typename ValueType>
void setCSRStorage( _MatrixStorage& storage )
{
    const IndexType numValues = 12;
    const IndexType numRows = 7;
    const IndexType numColumns = 4;
    ValueType values[] =
    { 6.0f, 4.0f, 7.0f, 9.0f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };
    IndexType ja[] =
    { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };
    IndexType ia[] =
    { 0, 2, 3, 5, 8, 10, 10, 12 };
    HArrayRef<IndexType> csrIA( numRows + 1, ia );
    HArrayRef<IndexType> csrJA( numValues, ja );
    HArrayRef<ValueType> csrValues( numValues, values );
    storage.setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
}

template<typename ValueType>
void setDenseStorage( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] =
    { 6.0f, 0.0f, 0.0f, 4.0f, 7.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 2.0f, 5.0f, 0.0f, 3.0f };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    storage.setRawDenseData( numRows, numColumns, values, ( ValueType ) 1E-5 );
}

template<typename ValueType1, typename ValueType2>
void checkEqual( MatrixStorage<ValueType1>& storage1, MatrixStorage<ValueType2>& storage2 )
{
    BOOST_REQUIRE_EQUAL( storage1.getNumRows(), storage2.getNumRows() );
    BOOST_REQUIRE_EQUAL( storage1.getNumColumns(), storage2.getNumColumns() );

    for ( IndexType i = 0; i < storage1.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < storage1.getNumColumns(); ++j )
        {
            float v1 = static_cast<float>( storage1.getValue( i, j ) );
            float v2 = static_cast<float>( storage2.getValue( i, j ) );
            BOOST_CHECK_CLOSE( v1, v2, 1.0e-5f );
        }
    }
}

template<typename StorageType1, typename StorageType2>
void conversion( ContextPtr loc )
{
    if ( !loc )
    {
        // skip this test
        return;
    }

    typedef typename StorageType1::StorageValueType ValueType1;
    // test empty assignment
    StorageType1 storage1;
    setCSRStorage<ValueType1>( storage1 );
    StorageType2 storage2( storage1, loc );
    checkEqual( storage1, storage2 );
    setDenseStorage<ValueType1>( storage1 );
    storage2 = storage1;
    checkEqual( storage1, storage2 );
}

/* ------------------------------------------------------------------------- */

static ContextPtr host = Context::getHostPtr();

static ContextPtr cuda =
    Context::hasContext( scai::common::context::CUDA ) ? Context::getContextPtr( scai::common::context::CUDA ) :
    ContextPtr();

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MatrixStorageConversionTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixStorageConversionTest" )

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CSR2ELLTest, ValueType, test_types )
{
    conversion< CSRStorage<ValueType>, ELLStorage<float> >( cuda );
    conversion< CSRStorage<ValueType>, ELLStorage<double> >( cuda );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( CSR2DenseTest, ValueType, test_types )
{
    conversion< CSRStorage<ValueType>, DenseStorage<float> >( host );
    conversion< CSRStorage<ValueType>, DenseStorage<double> >( host );
    conversion< CSRStorage<ValueType>, DenseStorage<float> >( cuda );
    conversion< CSRStorage<ValueType>, DenseStorage<double> >( cuda );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( Dense2CSRTest, ValueType, test_types )
{
    conversion< DenseStorage<ValueType>, CSRStorage<float> >( host );
    conversion< DenseStorage<ValueType>, CSRStorage<double> >( cuda );
}

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list< CSRStorage<double>,
        ELLStorage<float>,
        COOStorage<double>,
        JDSStorage<float>,
        DIAStorage<double>,
        DenseStorage<float>,
        SparseAssemblyStorage<double> > StorageTypes;

BOOST_AUTO_TEST_CASE_TEMPLATE( CreateTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    StorageType storage;
    setCSRStorage<ValueType>( storage );
    scai::common::shared_ptr<_MatrixStorage> storage1( storage.clone() );
// check for same format and value type
    BOOST_CHECK_EQUAL( storage1->getFormat(), storage.getFormat() );
    BOOST_CHECK_EQUAL( storage1->getValueType(), storage.getValueType() );
// check that it is a zero storage
    BOOST_CHECK_EQUAL( storage1->getNumRows(), 0 );
    BOOST_CHECK_EQUAL( storage1->getNumColumns(), 0 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( CopyTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    StorageType storage;
    setCSRStorage<ValueType>( storage );
    scai::common::shared_ptr<MatrixStorage<ValueType> > storage1( storage.copy() );
// check for same format and value type
    BOOST_CHECK_EQUAL( storage1->getFormat(), storage.getFormat() );
    BOOST_CHECK_EQUAL( storage1->getValueType(), storage.getValueType() );
// check that it has the same values
    checkEqual( storage, *storage1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( DefaultConstructorTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    StorageType storage;
    setCSRStorage<ValueType>( storage );
    StorageType storage1( storage );
// check for same format and value type
    BOOST_CHECK_EQUAL( storage1.getFormat(), storage.getFormat() );
    BOOST_CHECK_EQUAL( storage1.getValueType(), storage.getValueType() );
// check that it has the same values
    checkEqual( storage, storage1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    DenseStorage<ValueType> dense;
    setCSRStorage<ValueType>( dense );
    StorageType storage( dense );
    // check that it has the same values
    checkEqual( dense, storage );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( DefaultAssignmentTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    StorageType storage;
    setCSRStorage<ValueType>( storage );
    StorageType storage1;
    storage1 = storage;
// check that it has the same values
    checkEqual( storage, storage1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( AssignmentTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    CSRStorage<ValueType> csr;
    setCSRStorage<ValueType>( csr );
    StorageType storage;
    storage = csr;
// check that it has the same values
    checkEqual( csr, storage );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SwapTest, StorageType, StorageTypes )
{
    typedef typename StorageType::StorageValueType ValueType;
    StorageType storage, saveStorage;
    setCSRStorage<ValueType>( storage );
    setCSRStorage<ValueType>( saveStorage );
    StorageType storage1;
    storage1.swap( storage );
// check that storage is now a zero storage
    BOOST_CHECK_EQUAL( storage.getNumRows(), 0 );
    BOOST_CHECK_EQUAL( storage.getNumColumns(), 0 );
// check that storage1 now is the CSR matrix
    checkEqual( saveStorage , storage1 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

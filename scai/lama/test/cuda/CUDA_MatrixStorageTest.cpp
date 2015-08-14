/**
 * @file CUDAMatrixStorageTest.cpp
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
 * @brief Contains the implementation of the class CUDAMatrixStorageTest.cpp
 * @author: Alexander Büchel, Thomas Brandes, Bea Hornef
 * @date 02.05.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/lama/distribution/NoDistribution.hpp>

#include <test/cuda/CUDAContext.hpp>

using namespace lama;
using namespace memory;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_MatrixStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_CUDAMatrixStorageTest" );

/* --------------------------------------------------------------------- */

template<typename StorageType>
void setCSRData( StorageType& storage )
{
    typedef typename StorageType::StorageValueType ValueType; //!< This is the type of the matrix values.
    const IndexType numValues = 12;
    const IndexType numRows = 7;
    const IndexType numColumns = 4;
    static double values[] =
    { 6.0, 4.0, 7.0, 9.0, 4.0, 2.0, 5.0, 3.0, 2.0, 1.0, 1.0, 2.0 };
    static IndexType ja[] =
    { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };
    static IndexType ia[] =
    { 0, 2, 3, 5, 8, 10, 10, 12 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( int( sizeof( values ) / sizeof( double ) ), numValues );
    BOOST_CHECK_EQUAL( sizeof( values ) / sizeof( double ), sizeof( ja ) / sizeof( IndexType ) );
    BOOST_CHECK_EQUAL( int( sizeof( ia ) / sizeof( IndexType ) ), numRows + 1 );
    LAMAArray<IndexType> csrIas;
    LAMAArray<IndexType> csrJas;
    LAMAArray<ValueType> csrValues;
    WriteOnlyAccess<IndexType> myIa( csrIas, numRows + 1 );
    WriteOnlyAccess<IndexType> myJa( csrJas, numValues );
    WriteOnlyAccess<ValueType> myData( csrValues, numValues );

    for ( IndexType ii = 0; ii <= numRows; ii++ )
    {
        myIa[ii] = ia[ii];
    }

    for ( IndexType jj = 0; jj < numValues; jj++ )
    {
        myJa[jj] = ja[jj];
        myData[jj] = static_cast<ValueType>( values[jj] );
    }

    myIa.release();
    myJa.release();
    myData.release();
    ContextPtr host = Context::getContextPtr( context::Host );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    storage.setContext( host );
    storage.setCSRData( numRows, numColumns, numValues, csrIas, csrJas, csrValues );

    // fill with the csr sparse data
    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            IndexType j = ja[jj];
            BOOST_CHECK_CLOSE( values[jj], storage.getValue( i, j ), 1e-16 );
        }
    }

    storage.setContext( cuda );
    storage.setCSRData( numRows, numColumns, numValues, csrIas, csrJas, csrValues );

    // fill with the csr sparse data
    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            IndexType j = ja[jj];
            BOOST_CHECK_CLOSE( values[jj], storage.getValue( i, j ), 1e-16 );
        }
    }
}

template<typename StorageType>
void testSetMethod()
{
    StorageType storage;
    setCSRData( storage );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( setXToCSRData, ValueType, test_types )
{
    testSetMethod<ELLStorage<ValueType> >();
    testSetMethod<DenseStorage<ValueType> >();
    testSetMethod<CSRStorage<ValueType> >();
    testSetMethod<COOStorage<ValueType> >();
    testSetMethod<DIAStorage<ValueType> >();
    testSetMethod<JDSStorage<ValueType> >();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

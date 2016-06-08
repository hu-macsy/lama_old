/**
 * @file CUDA_MatrixStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Contains the implementation of the class CUDAMatrixStorageTest.cpp
 * @author Alexander Büchel, Thomas Brandes, Bea Hornef
 * @date 02.05.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/dmemo/NoDistribution.hpp>

#include <scai/lama/test/cuda/CUDATestContext.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

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
    HArray<IndexType> csrIas;
    HArray<IndexType> csrJas;
    HArray<ValueType> csrValues;
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
    ContextPtr host = Context::getHostPtr();
    ContextPtr cuda = scai::lama_test::CUDATestContext::getContext();
    storage.setContextPtr( host );
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

    storage.setContextPtr( cuda );
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

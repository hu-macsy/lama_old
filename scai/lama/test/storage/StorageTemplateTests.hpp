/**
 * @file StorageTemplateTests.hpp
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
 * @brief Storage tests that require the storage type as template param
 * @author Thomas Brandes
 * @date 20.04.2015
 */

#pragma once

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

/** Test of swap method that is defined only for two storages of same format, type. */

template<typename StorageType>
void storageSwapTest()
{
    using namespace scai;
    using namespace hmemo;
    using namespace utilskernel;
    using namespace lama;

    typedef typename StorageType::StorageValueType ValueType;

    ContextPtr context = Context::getContextPtr();

    StorageType csr1;
    StorageType csr2;

    setDenseData( csr1 );

    IndexType n = csr1.getNumRows();
    IndexType m = csr1.getNumColumns();

    LArray<ValueType> x( context );
    LArray<ValueType> y( context );
    LArray<ValueType> z1( context );
    LArray<ValueType> z2( context );

    ValueType alpha = 1.3;
    ValueType beta  = -0.5;

    HArrayUtils::setRandom( x, m, 1.0f );
    HArrayUtils::setRandom( y, n, 1.0f );

    csr1.matrixTimesVector( z1, alpha, x, beta, y );

    csr1.swap( csr2 );

    // now check sizes 

    BOOST_CHECK_EQUAL( n, csr2.getNumRows() );
    BOOST_CHECK_EQUAL( m, csr2.getNumColumns() );

    BOOST_CHECK_EQUAL( 0, csr1.getNumRows() );
    BOOST_CHECK_EQUAL( 0, csr1.getNumColumns() );

    // now check that the other matrix contains the right values

    csr2.matrixTimesVector( z2, alpha, x, beta, y );

    typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

    AbsType diff = common::Math::real( z1.maxDiffNorm( z2 ) );

    // even if storages are the same we can have different rounding errors

    BOOST_CHECK( diff < 0.001 );
}

/** Each storage class has a static method typeName that is tested here. 
 *
 *  @param[in] subString  this string must be in the typename
 */

template<typename StorageType>
void storageTypeNameTest( const char* subString )
{
    // context does not matter here, so runs for every context

    std::string s = StorageType::typeName();

    BOOST_CHECK( s.length() > 0 );

    // substring must be in typename

    BOOST_CHECK( s.find( subString ) != std::string::npos );
}

template<typename StorageType>
void copyStorageTest()
{
    StorageType storage;

    typedef typename StorageType::StorageValueType ValueType;

    float defaultThreshold = storage.getCompressThreshold();
    float fullThreshold = 1.0;

    BOOST_CHECK( fullThreshold != defaultThreshold );

    storage.setCompressThreshold( fullThreshold );

    setDenseHalo( storage );

    const scai::lama::MatrixStorage<ValueType>& matrixStorage = storage;

    StorageType copyStorage1( storage );          // default copy constructor
    StorageType copyStorage2( matrixStorage );    // own copy constructor
    StorageType copyStorage3;                     // default assignment operator
    copyStorage3 = storage;
    StorageType copyStorage4;                     // own assignment operator
    copyStorage4 = matrixStorage;

    // if default constructors are not overwritten compressThreshold is also copied

    BOOST_CHECK_EQUAL( copyStorage1.getCompressThreshold(), defaultThreshold );
    BOOST_CHECK_EQUAL( copyStorage2.getCompressThreshold(), defaultThreshold );
    BOOST_CHECK_EQUAL( copyStorage3.getCompressThreshold(), defaultThreshold );
    BOOST_CHECK_EQUAL( copyStorage4.getCompressThreshold(), defaultThreshold );
}


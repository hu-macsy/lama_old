/**
 * @file StorageTemplateTests.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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

namespace scai
{

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
    const lama::MatrixStorage<ValueType>& matrixStorage = storage;

    StorageType copyStorage1( storage );          // default copy constructor
    StorageType copyStorage2 = lama::convert<StorageType>( matrixStorage );    // own copy constructor
    StorageType copyStorage3;                     // default assignment operator
    copyStorage3 = storage;
    StorageType copyStorage4;                     // own assignment operator
    copyStorage4.assign( matrixStorage );

    // copy constructors take over the threshold

    BOOST_CHECK_EQUAL( copyStorage1.getCompressThreshold(), fullThreshold );

    // assignment operator do not modify the threshold

    BOOST_CHECK_EQUAL( copyStorage2.getCompressThreshold(), defaultThreshold );
    BOOST_CHECK_EQUAL( copyStorage3.getCompressThreshold(), defaultThreshold );
    BOOST_CHECK_EQUAL( copyStorage4.getCompressThreshold(), defaultThreshold );
}

}

/**
 * @file StorageIOTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class StorageIOTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/StorageIO.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai::common;
using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StorageIOTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.StorageIOTest" );

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 7;
    const IndexType numColumns = 7;
    static ValueType values[] =
    {
        6.0f, 0.0f, 0.0f, 4.0f, 0.0f, 1.0f, -2.11f, 7.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 3.0f, 1.0f, -1.0f, 0.0f, 0.0f, 1.1f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 2.3f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 1.0f, -3.0f
    };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOFormatted, ValueType, scai_arithmetic_test_types )
{
    if ( isComplex( TypeTraits<ValueType>::stype ) )
    {
        return;
    }

    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_formatted";
    csrMatrix.writeToFile( filename, File::SAMG_FORMAT );
    readMatrix.readFromFile( filename + ".frm" );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumRows(), csrMatrix.getNumRows() );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumColumns(), csrMatrix.getNumColumns() );

    for ( IndexType i = 0; i < csrMatrix.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrMatrix.getNumColumns(); ++j )
        {
            SCAI_CHECK_CLOSE( csrMatrix.getValue( i, j ),
                              readMatrix.getValue( i, j ), 0.01f );
        }
    }

    std::remove( ( filename + ".amg" ).c_str() );
    std::remove( ( filename + ".frm" ).c_str() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOmatrixMarket, ValueType, scai_arithmetic_test_types )
{
    scalar::ScalarType stype = TypeTraits<ValueType>::stype;

    if ( isComplex( stype ) || ( stype == scalar::LONG_DOUBLE ) )
    {
        return;
    }

    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_matrix_market";
    csrMatrix.writeToFile( filename, File::MATRIX_MARKET );
    readMatrix.readFromFile( filename + ".mtx" );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumRows(), csrMatrix.getNumRows() );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumColumns(), csrMatrix.getNumColumns() );

    for ( IndexType i = 0; i < csrMatrix.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrMatrix.getNumColumns(); ++j )
        {
            SCAI_CHECK_CLOSE( csrMatrix.getValue( i, j ), readMatrix.getValue( i, j ), 0.01f );
        }
    }

    std::remove( filename.c_str() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOBinary, ValueType, scai_arithmetic_test_types )
{
    if ( isComplex( TypeTraits<ValueType>::stype ) )
    {
        return;
    }

    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_binary";
    csrMatrix.writeToFile( filename, File::SAMG_FORMAT, TypeTraits<ValueType>::stype, TypeTraits<IndexType>::stype, TypeTraits<IndexType>::stype, true );
    readMatrix.readFromFile( filename + ".frm" );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumRows(), csrMatrix.getNumRows() );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumColumns(), csrMatrix.getNumColumns() );

// Binary data has no lost of accuracy, so we compare for equality here

    for ( IndexType i = 0; i < csrMatrix.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrMatrix.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( csrMatrix.getValue( i, j ),
                               readMatrix.getValue( i, j ) );
        }
    }

    // TODO: uncomment
//    std::remove( ( filename + ".amg" ).c_str() );
//    std::remove( ( filename + ".frm" ).c_str() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

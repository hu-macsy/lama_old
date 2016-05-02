/**
 * @file StorageIOTest.cpp
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
 * @brief Contains the implementation of the class StorageIOTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/StorageIO.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Complex.hpp>

using namespace scai::common;
using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StorageIOTest )
;

SCAI_LOG_DEF_LOGGER( logger, "Test.StorageIOTest" );

typedef boost::mpl::list<float, double> test_types;
typedef boost::mpl::list<float, double, Complex<float>, Complex<double> > test_types_complex;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOFormatted, ValueType, test_types )
{
    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_formatted";
    csrMatrix.writeToFile( filename, File::SAMG );
    readMatrix.readFromFile( filename + ".frm" );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumRows(), csrMatrix.getNumRows() );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumColumns(), csrMatrix.getNumColumns() );

    for ( IndexType i = 0; i < csrMatrix.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrMatrix.getNumColumns(); ++j )
        {
            BOOST_CHECK_CLOSE( csrMatrix.getValue( i, j ),
                               readMatrix.getValue( i, j ), 0.01f );
        }
    }

    std::remove( ( filename + ".amg" ).c_str() );
    std::remove( ( filename + ".frm" ).c_str() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOmatrixMarket, ValueType, test_types_complex )
{
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

BOOST_AUTO_TEST_CASE_TEMPLATE( StorageIOBinary, ValueType, test_types )
{
    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_binary";
    csrMatrix.writeToFile( filename, File::SAMG, TypeTraits<ValueType>::stype, TypeTraits<IndexType>::stype, TypeTraits<IndexType>::stype, true );
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

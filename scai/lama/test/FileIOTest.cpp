/**
 * @file FileIOTest.cpp
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
 * @brief Contains the implementation of the class FileIOTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace scai::common;
using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( FileIOTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.FileIOTest" );

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 7;
    const IndexType numColumns = 7;

    // values: take numRows x numColums random numbers of required type

    HArray<ValueType> values;
    float fillRate = 0.2;    
    utilskernel::HArrayUtils::setRandom( values, numRows * numColumns, fillRate );

    ValueType eps = static_cast<ValueType>( 1E-5 );

    // Note: diagonal property of sparse matrices will be set due to square matrix

    storage.setDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FileIOFormatted, ValueType, scai_arithmetic_test_types )
{
    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
 
    SCAI_LOG_INFO( logger, "FileIOFormatted: write SAMG Formatted of this storage: " << csrMatrix )

    std::string filename = "out_formatted.frm";
    csrMatrix.writeToFile( filename, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::FORMATTED );
    BOOST_CHECK( FileIO::fileExists( filename ) );
    readMatrix.readFromFile( filename );
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

    int rc = FileIO::removeFile( filename );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( filename ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FileIOMatrixMarket, ValueType, scai_arithmetic_test_types )
{
    // scalar::ScalarType stype = TypeTraits<ValueType>::stype;

    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );

    std::string filename = "out_matrix_market.mtx";   

    SCAI_LOG_INFO( logger, "FileIOMatrixMarket: write to file " << filename << " this storage: " << csrMatrix )

    // specify all values to ignore any settings
    csrMatrix.writeToFile( filename, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::FORMATTED );

    BOOST_REQUIRE( FileIO::fileExists( filename ) );

    readMatrix.readFromFile( filename );

    BOOST_REQUIRE_EQUAL( readMatrix.getNumRows(), csrMatrix.getNumRows() );
    BOOST_REQUIRE_EQUAL( readMatrix.getNumColumns(), csrMatrix.getNumColumns() );

    for ( IndexType i = 0; i < csrMatrix.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrMatrix.getNumColumns(); ++j )
        {
            SCAI_CHECK_CLOSE( csrMatrix.getValue( i, j ), readMatrix.getValue( i, j ), 0.01f );
        }
    }

    int rc = FileIO::removeFile( filename );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( filename ) );

    BOOST_CHECK_THROW (
    {
        readMatrix.readFromFile( filename );
    },
    Exception );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FileIOBinary, ValueType, scai_arithmetic_test_types )
{
    CSRStorage<ValueType> csrMatrix;
    CSRStorage<ValueType> readMatrix;
    setDenseData( csrMatrix );
    std::string filename = "out_binary.frm";

    csrMatrix.writeToFile( filename, "", scalar::INTERNAL, scalar::INDEX_TYPE, FileIO::BINARY );

    BOOST_REQUIRE( FileIO::fileExists( filename ) );

    readMatrix.readFromFile( filename );
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

    int rc = FileIO::removeFile( filename );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( filename ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

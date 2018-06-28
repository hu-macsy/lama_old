/**
 * @file MMIOTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Individual tests for the matrix market format
 * @author Thomas Brandes
 * @date 03.01.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/io/MatrixMarketIO.hpp>
#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/macros/assert.hpp>

#include <scai/testsupport/uniquePath.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace common;
using namespace lama;
using namespace hmemo;

using scai::testsupport::uniquePath;
using scai::testsupport::GlobalTempDir;

/** Output files should be deleted unless for debugging it might be useful to check them. */

#undef DELETE_OUTPUT_FILES

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MMIOTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MMIOTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ReadGeneralDenseTest, ValueType, scai_numeric_test_types )
{
    const IndexType m   = 3;
    const IndexType n   = 4;

    ValueType vals[]  = { 1, 2, 3, 0, 1, -1, 2, 0, 1, 5, 2, 7 };

    const IndexType n_vals = sizeof( vals ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( n_vals, m * n, "number of values must fit to matrix " << m << " x " << n );

    const char header[] = "%%MatrixMarket matrix array real general";

    const auto fileName = uniquePath(GlobalTempDir::getPath(), "mm_array_gen") + ".mtx";

    {
        using namespace std;

        fstream myFile;

        myFile.open( fileName.c_str(), ios::out );
        myFile << header << endl;
        myFile << m << " " << n << endl;

        for ( IndexType i = 0; i < m * n; ++i )
        {
            myFile << vals[i] << endl;
        }
    }

    DenseStorage<ValueType> denseStorage;

    MatrixMarketIO reader;
    FileIO& freader = reader;

    freader.readStorage( denseStorage, fileName );

    BOOST_CHECK_EQUAL( m, denseStorage.getNumRows() );
    BOOST_CHECK_EQUAL( n, denseStorage.getNumColumns() );

    {
        const HArray<ValueType>& data = denseStorage.getData();
        ReadAccess<ValueType> rVals( data );

        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = 0; j < n; ++j )
            {
                // Note: denseStorage stores row-wise, file was column-wise
                // std::cout << "Array [ " << i << ", " << j << " ] = " << rVals[ i * n + j ] << std::endl;

                BOOST_CHECK_EQUAL( rVals[i * n + j ], vals[ j * m + i ] );
            }
        }
    }

    int rc = FileIO::removeFile( fileName );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( fileName ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ReadErrorTest )
{
    using namespace std;

    const char header_wrong1[] = "%%MatMarket matrix array real general";
    const char header_wrong2[] = "%%MatrixMarket mat array real general";
    const char header_wrong3[] = "%%MatrixMarket vector coordinat real symmetric";
    const char header_wrong4[] = "%%MatrixMarket matrix coordinate doubl general";
    const char header_wrong5[] = "%%MatrixMarket matrix coordinate real hermitan";

    const char* header[] = { header_wrong1, header_wrong2, header_wrong3, header_wrong4, header_wrong5 };

    const IndexType ncases = sizeof( header ) / sizeof( char* );

    for ( IndexType icase = 0; icase < ncases; ++icase )
    {
        const auto fileName = uniquePath(GlobalTempDir::getPath(), "mm_wrong") + ".mtx";

        {
            fstream myFile;
            IndexType zero = 0;    // used to write sizes

            myFile.open( fileName.c_str(), ios::out );
            myFile << header[icase] << endl;
            myFile << zero << " " << zero << endl;
        }

        MatrixMarketIO reader;
 
        SCAI_LOG_INFO( logger, "Test this header: " << header[icase] )

        BOOST_CHECK_THROW(
        {
            IndexType N = 0;
            reader.readArrayInfo( N, fileName );
        }, common::IOException );

        int rc = FileIO::removeFile( fileName );

        BOOST_CHECK_EQUAL( rc, 0 );
        BOOST_CHECK( ! FileIO::fileExists( fileName ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ReadSymmetricDenseTest, ValueType, scai_numeric_test_types )
{
    const IndexType n   = 4;

    /*    DenseMatrix    1  -  -  -
                         2  1  -  -
                         3 -1  5  -
                         0 -2  7  6
    */

    ValueType vals[]  = { 1, 2, 3, 0, 1, -1, -2, 5, 7, 6 };

    const IndexType n_vals = sizeof( vals ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( n_vals, n * ( n + 1 ) / 2 , "number of symmentric values must fit to matrix " << n << " x " << n );

    const char header[] = "%%MatrixMarket matrix array real symmetric";

    const auto fileName = uniquePath(GlobalTempDir::getPath(), "mm_array_symm") + ".mtx";

    {
        using namespace std;

        fstream myFile;

        myFile.open( fileName.c_str(), ios::out );
        myFile << header << endl;
        myFile << n << " " << n << endl;

        for ( IndexType i = 0; i < n_vals; ++i )
        {
            myFile << vals[i] << endl;
        }
    }

    DenseStorage<ValueType> denseStorage;

    MatrixMarketIO reader;
    FileIO& freader = reader;

    freader.readStorage( denseStorage, fileName );

    BOOST_CHECK_EQUAL( n, denseStorage.getNumRows() );
    BOOST_CHECK_EQUAL( n, denseStorage.getNumColumns() );

    {
        const HArray<ValueType>& data = denseStorage.getData();
        ReadAccess<ValueType> rVals( data );

        IndexType k = 0;

        for ( IndexType i = 0; i < n; ++i )
        {
            for ( IndexType j = i; j < n; ++j )
            {
                //  upper and lower must be same

                BOOST_CHECK_EQUAL( rVals[i * n + j ], vals[ k ] );
                BOOST_CHECK_EQUAL( rVals[j * n + i ], vals[ k ] );

                k++;
            }
        }

        SCAI_ASSERT_EQ_ERROR( k, n_vals, "mismatch for number of values traversing triangular matrix" )
    }

    int rc = FileIO::removeFile( fileName );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( fileName ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( RectangularTest )
{
    const IndexType m   = 3;
    const IndexType n   = 5;

    DefaultReal vals[]  = { 1, 2, 3,
                            4, 5, 6,
                            7, 8, 9,
                            10, 11, 12,
                            13, 14, 15 };

    DenseStorage<DefaultReal> denseOut;
    DenseStorage<DefaultReal> denseIn;

    denseOut.setRawDenseData( m, n, vals );

    const auto fileName = uniquePath(GlobalTempDir::getPath(), "dense") + ".mtx";

    MatrixMarketIO io;
    FileIO& fio = io;

    fio.writeStorage( denseOut, fileName );
    fio.readStorage( denseIn, fileName );

    BOOST_CHECK_EQUAL( denseOut.getNumColumns(), denseIn.getNumColumns() );
    BOOST_CHECK_EQUAL( denseIn.getNumRows(), denseIn.getNumRows() );

    int rc = FileIO::removeFile( fileName );

    BOOST_CHECK_EQUAL( rc, 0 );
    BOOST_CHECK( ! FileIO::fileExists( fileName ) );

    DefaultReal diff = denseIn.maxDiffNorm( denseOut );
    BOOST_CHECK( diff < 0.0001 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */

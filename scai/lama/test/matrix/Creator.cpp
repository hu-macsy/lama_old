/**
 * @file lama/test/matrix/Creator.cpp
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
 * @brief Create distributed matrices
 * @author Thomas Brandes
 * @date 31.05.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace dmemo;

/* ------------------------------------------------------------------------- */

// This routine fills a matrix storage with some data

template<typename ValueType>
static void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 2;
    const IndexType numColumns = 8;
    static ValueType values[] =
    { 6.0, 0.0, 0.0, 4.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 4.0, 2.0, 5.0, 0.0, 3.0 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( Creator )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.Create" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildReplicatedDiagTest )
{
    // it should be okay just to test one ValueType

    typedef SCAI_TEST_TYPE ValueType;

    CSRStorage<ValueType> storage;
    setDenseData( storage );

    SCAI_LOG_DEBUG( logger, "storage to replicate: " << storage )

    CSRSparseMatrix<ValueType> matrix;

    const IndexType repN = 6;

    MatrixCreator::buildReplicatedDiag( matrix, storage, repN );

    SCAI_LOG_INFO( logger, matrix.getRowDistribution().getReduceCommunicator() << ": buildReplicatedDiag: " << matrix
                   << " from rep = " << repN << " x " << storage )

    const IndexType nRows = storage.getNumRows();
    const IndexType nCols = storage.getNumColumns();

    // first check for correct sizes

    BOOST_CHECK_EQUAL( matrix.getNumRows(), repN * nRows );
    BOOST_CHECK_EQUAL( matrix.getNumColumns(), repN * nCols );

    // traverse matrix row-wise in blocks to check correct values

    for ( IndexType i1 = 0; i1 < repN; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < nRows; ++i2 )
        {
            IndexType i = i1 * nRows + i2;

            for ( IndexType j1 = 0; j1 < repN; ++j1 )
            {
                for ( IndexType j2 = 0; j2 < nCols; ++j2 )
                {
                    IndexType j = j1 * nCols + j2;

                    ValueType v = matrix.getValue( i, j );

                    SCAI_LOG_TRACE( logger, "i = " << i << " ( " << i1 << " * " << nRows << " + " << i2 << " )"
                                    << ", j = " << j << " ( " << j1 << " * " << nCols << " + " << i2 << " )"
                                    << ", v = " << v )

                    if ( i1 != j1 )
                    {
                        // not a diagonal block, so it must be all 0

                        BOOST_CHECK_EQUAL( v, 0.0 );
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( v, storage.getValue( i2, j2 ) );
                    }
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildReplicatedTest )
{
    // it should be okay just to test one ValueType

    typedef SCAI_TEST_TYPE ValueType;

    CSRStorage<ValueType> storage;
    setDenseData( storage );

    const IndexType nRows = storage.getNumRows();
    const IndexType nCols = storage.getNumColumns();

    CSRSparseMatrix<ValueType> matrix;

    const IndexType repRow = 2;
    const IndexType repCol = 2;

    MatrixCreator::buildReplicated( matrix, storage, repRow, repCol );

    SCAI_LOG_INFO( logger, matrix.getRowDistribution().getReduceCommunicator() << ": buildReplicated " << matrix )

    // first check for correct sizes

    BOOST_CHECK_EQUAL( matrix.getNumRows(), repRow * nRows );
    BOOST_CHECK_EQUAL( matrix.getNumColumns(), repCol * nCols );

    // traverse matrix row-wise in blocks to check correct values

    for ( IndexType i1 = 0; i1 < repRow; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < nRows; ++i2 )
        {
            IndexType i = i1 * nRows + i2;

            for ( IndexType j1 = 0; j1 < repCol; ++j1 )
            {
                for ( IndexType j2 = 0; j2 < nCols; ++j2 )
                {
                    IndexType j = j1 * nCols + j2;

                    ValueType v = matrix.getValue( i, j );
                    BOOST_CHECK_EQUAL( v, storage.getValue( i2, j2 ) );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */

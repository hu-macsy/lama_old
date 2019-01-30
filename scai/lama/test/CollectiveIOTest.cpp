/**
 * @file CollectioveIOTest.cpp
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
 * @brief Contains specific tests for class CollectioveIO 
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/CollectiveIO.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/testsupport/uniquePathComm.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace lama;

using scai::testsupport::uniquePathSharedAmongNodes;
using scai::testsupport::GlobalTempDir;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CollectiveIOTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CollectiveIOTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DenseVectorTest, ValueType, scai_array_test_types )
{
    using namespace hmemo;

    const IndexType N = 30;

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    auto dist = dmemo::blockDistribution( N, comm );

    ValueType inc = ValueType( 1 ) / ValueType( N );

    auto v1 = denseVectorLinear<ValueType>( dist, 0, inc );

    const char fileName[] = "DenseVectorTest.data";

    auto file = comm->collectiveFile();

    file->open( fileName, "w" );
    CollectiveIO::write( *file, v1 );
    file->close();
    
    file->open( fileName, "r" );
    DenseVector<ValueType> v2;
    CollectiveIO::read( *file, v2 );
    file->close();

    BOOST_TEST( hostReadAccess( v1.getLocalValues() ) == hostReadAccess( v2.getLocalValues() ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName );
    };
}

/* --------------------------------------------------------------------- */

/**
 *  Help function to write arbirtrary dense vectors in a collective file.
 */
template<typename FileIndexType, typename FileDataType>
static void writeTestDenseVector( dmemo::CollectiveFile& outFile, const IndexType N )
{
    auto indexType = common::TypeTraits<FileIndexType>::stype;
    auto dataType = common::TypeTraits<FileDataType>::stype;

    outFile.writeSingle( CollectiveIO::getDenseVectorId() );
    outFile.writeSingle( static_cast<int>( indexType ) );
    outFile.writeSingle( static_cast<int>( dataType ) );
    outFile.writeSingle( N, indexType );

    for ( IndexType i = 0; i < N; ++i )
    {
        outFile.writeSingle( i, dataType );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DenseVectorConvertTest, ValueType, scai_array_test_types )
{
    using namespace hmemo;

    const IndexType N = 10;

    auto comm = dmemo::Communicator::getCommunicatorPtr();
    auto dist = dmemo::blockDistribution( N, comm );

    const std::string fileName = uniquePathSharedAmongNodes(
                                      GlobalTempDir::getPath(),
                                      *comm,
                                      "Test"
                                    ) + ".raw";

    auto file = comm->collectiveFile();

    file->open( fileName.c_str(), "w" );
    writeTestDenseVector<int, DefaultReal>( *file, N );
    writeTestDenseVector<unsigned long, IndexType>( *file, N );
    writeTestDenseVector<float, DefaultReal>( *file, N );
    file->close();
    
    DenseVector<ValueType> v1, v2, v3;

    file->open( fileName.c_str(), "r" );
    CollectiveIO::read( *file, v1 );
    CollectiveIO::read( *file, v2 );
    CollectiveIO::read( *file, v3 );
    file->close();

    auto expV = denseVectorLinear<ValueType>( dist, 0, 1 );

    BOOST_TEST( hostReadAccess( expV.getLocalValues() ) == hostReadAccess( v1.getLocalValues() ), per_element() );
    BOOST_TEST( hostReadAccess( expV.getLocalValues() ) == hostReadAccess( v2.getLocalValues() ), per_element() );
    BOOST_TEST( hostReadAccess( expV.getLocalValues() ) == hostReadAccess( v3.getLocalValues() ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName.c_str() );
    };
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SparseVectorTest )
{
    typedef DefaultReal ValueType;

    using namespace hmemo;

    const IndexType N = 30;

    HArray<IndexType> indexes( { 0, 5, 17, 21 } );
    HArray<ValueType> values( { 1, 2, 3, -1 } );

    auto comm = dmemo::Communicator::getCommunicatorPtr();
    auto dist = dmemo::blockDistribution( N, comm );

    SparseVector<ValueType> testSV( N, indexes, values );
    testSV.redistribute( blockDistribution( N, comm ) );

    const char fileName[] = "SparseVectorTest.RAW";

    auto file = comm->collectiveFile();

    file->open( fileName, "w" );
    CollectiveIO::write( *file, testSV );
    file->close();
    
    SparseVector<ValueType> readSV;

    file->open( fileName, "r" );
    CollectiveIO::read( *file, readSV );
    file->close();

    SCAI_LOG_INFO( logger, "read sparse vector from file " << fileName << ": " << readSV )

    readSV.redistribute( testSV.getDistributionPtr() );

    BOOST_TEST( hostReadAccess( testSV.getNonZeroIndexes() ) == hostReadAccess( readSV.getNonZeroIndexes() ), per_element() );
    BOOST_TEST( hostReadAccess( testSV.getNonZeroValues() ) == hostReadAccess( readSV.getNonZeroValues() ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName );
    };
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CSRSparseMatrixTest )
{
    using namespace hmemo;

    typedef DefaultReal ValueType;

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType numRows = 15;   // # rows for global matrix
    const IndexType numCols = 10;   // # cols for global matrix

    auto rowDist = dmemo::blockDistribution( numRows );
    auto colDist = dmemo::noDistribution( numCols );

    const std::string fileName  = "csrMatrix.RAW";

    auto matrix = zero<CSRSparseMatrix<ValueType>>( rowDist, colDist );

    float fillRate = 0.2f;
    MatrixCreator::fillRandom( matrix, fillRate );

    auto file = comm->collectiveFile();

    SCAI_LOG_INFO( logger, "write CSR sparse matrix into file " << fileName );

    file->open( fileName.c_str(), "w" );
    CollectiveIO::write( *file, matrix );
    file->close();
    
    CSRSparseMatrix<ValueType> readMatrix;

    SCAI_LOG_INFO( logger, "read CSR sparse matrix from file " << fileName );

    file->open( fileName.c_str(), "r" );
    CollectiveIO::read( *file, readMatrix );
    file->close();

    const auto& expIA = matrix.getLocalStorage().getIA();
    const auto& csrIA = readMatrix.getLocalStorage().getIA();

    const auto& expJA = matrix.getLocalStorage().getJA();
    const auto& csrJA = readMatrix.getLocalStorage().getJA();

    const auto& expValues = matrix.getLocalStorage().getValues();
    const auto& csrValues = readMatrix.getLocalStorage().getValues();

    BOOST_TEST( hostReadAccess( csrIA ) == hostReadAccess( expIA ), per_element() );
    BOOST_TEST( hostReadAccess( csrJA ) == hostReadAccess( expJA ), per_element() );
    BOOST_TEST( hostReadAccess( csrValues ) == hostReadAccess( expValues ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName.c_str() );
    };
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DenseMatrixTest )
{
    using namespace hmemo;

    typedef DefaultReal ValueType;

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType numRows = 15;   // # rows for global matrix
    const IndexType numCols = 10;   // # cols for global matrix

    auto rowDist = dmemo::blockDistribution( numRows );
    auto colDist = dmemo::noDistribution( numCols );

    const std::string fileName  = "denseMatrix.RAW";

    auto matrix = zero<DenseMatrix<ValueType>>( rowDist, colDist );

    float fillRate = 1.0f;
    MatrixCreator::fillRandom( matrix, fillRate );

    auto file = comm->collectiveFile();

    SCAI_LOG_INFO( logger, "write dense matrix into file " << fileName );

    file->open( fileName.c_str(), "w" );
    CollectiveIO::write( *file, matrix );
    file->close();
    
    DenseMatrix<ValueType> readMatrix;

    SCAI_LOG_INFO( logger, "read dense matrix from file " << fileName );

    file->open( fileName.c_str(), "r" );
    CollectiveIO::read( *file, readMatrix );
    file->close();

    const auto& expValues = matrix.getLocalStorage().getValues();
    const auto& denseValues = readMatrix.getLocalStorage().getValues();

    BOOST_TEST( hostReadAccess( denseValues ) == hostReadAccess( expValues ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName.c_str() );
    };
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

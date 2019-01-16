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

    auto v1 = linearDenseVector<ValueType>( dist, 0, inc );

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
    outFile.writeSingle( CollectiveIO::getDenseVectorId() );
    outFile.writeSingle( static_cast<int>( common::TypeTraits<FileIndexType>::stype ) );
    outFile.writeSingle( static_cast<int>( common::TypeTraits<FileDataType>::stype ) );
    outFile.writeSingle( static_cast<FileIndexType>( N ) );

    for ( IndexType i = 0; i < N; ++i )
    {
        outFile.writeSingle( static_cast<FileDataType>( i ) );
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
    BOOST_CHECK_THROW(
        {
            CollectiveIO::read( *file, v3 );  // float is no supported IndexType
        }, common::Exception );
    file->close();

    auto expV = linearDenseVector<ValueType>( dist, 0, 1 );

    BOOST_TEST( hostReadAccess( expV.getLocalValues() ) == hostReadAccess( v1.getLocalValues() ), per_element() );
    BOOST_TEST( hostReadAccess( expV.getLocalValues() ) == hostReadAccess( v2.getLocalValues() ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName.c_str() );
    };
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();


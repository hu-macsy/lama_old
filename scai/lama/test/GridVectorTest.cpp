/**
 * @file GridVectorTest.cpp
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
 * @brief Contains specific tests for class GridVector 
 * @author Thomas Brandes
 * @date 16.05.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/testsupport/uniquePathComm.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace lama;

using scai::testsupport::uniquePathSharedAmongNodes;
using scai::testsupport::GlobalTempDir;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GridVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GridVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_numeric_test_types )
{
    // replicated grid vector of a certain size only with zero elements

    const IndexType n1 = 4;
    const IndexType n2 = 5;
    const IndexType n3 = 6;

    const common::Grid3D grid( n1, n2, n3 );

    GridVector<ValueType> gv( grid, ValueType( 1 ) );

    BOOST_CHECK_EQUAL( 3, gv.nDims() );
    BOOST_CHECK_EQUAL( n1, gv.size( 0 ) );
    BOOST_CHECK_EQUAL( n2, gv.size( 1 ) );
    BOOST_CHECK_EQUAL( n3, gv.size( 2 ) );

    // be careful: gv.size( 3 ) only throws exception with SCAI_ASSERT_DEBUG enabled:x
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ioTest, ValueType, scai_numeric_test_types )
{
    if ( !FileIO::canCreate( ".mat" ) )
    {
        // Matlab IO might be disabled if no zlib was available
        return;
    }

    const IndexType n1 = 5;
    const IndexType n2 = 3;
    const IndexType n3 = 2;
    const IndexType n4 = 4;

    // Matlab binary format supports io of grids
    
    const auto comm = scai::dmemo::Communicator::getCommunicatorPtr();
    const auto fileName = uniquePathSharedAmongNodes(GlobalTempDir::getPath(), *comm, "tmpGrid") + ".mat";

    GridVector<ValueType> gv1( common::Grid4D( n1, n2, n3, n4 ), ValueType( 0 ) );

    {
        GridWriteAccess<ValueType> wGV1( gv1 );

        for ( IndexType i1 = 0; i1 < n1; ++i1 )
        {
            for ( IndexType i2 = 0; i2 < n2; ++i2 )
            {
                for ( IndexType i3 = 0; i3 < n3; ++i3 )
                {
                    for ( IndexType i4 = 0; i4 < n4; ++i4 )
                    { 
                        wGV1( i1, i2, i3, i4 ) = 
                            static_cast<ValueType>( 1000 * ( i1 + 1 ) + 100 * ( i2 + 1 ) + 10 * ( i3 + 1 ) + i4 + 1 );
                    }
                }
            }
        }
    }

    // write the grid vector to a file, might be done by first processor only

    gv1.writeToFile( fileName );

    // read the grid vector from a file

    auto gv2 = read<GridVector<ValueType>>( fileName );

    SCAI_LOG_INFO( logger, "gv2( " << fileName << " ) = " << gv2 )

    BOOST_CHECK_EQUAL( gv1.globalGrid(), gv2.globalGrid() );

    {
        GridReadAccess<ValueType> rGV1( gv1 );
        GridReadAccess<ValueType> rGV2( gv2 );

        for ( IndexType i = 0; i < gv1.size(); ++i )
        {
            SCAI_ASSERT_EQ_ERROR( rGV1[i], rGV2[i], "different val at i = " << i )
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

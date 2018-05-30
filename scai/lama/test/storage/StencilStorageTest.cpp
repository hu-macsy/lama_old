/**
 * @file StencilStorageTest.cpp
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
 * @brief Contains some test of the template class StencilStorage
 * @author Thomas Brandes
 * @date 24.05.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/StencilStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/utilskernel.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

using boost::test_tools::per_element;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.StencilStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    const IndexType N1 = 5;
    const IndexType N2 = 8;
    common::Grid2D grid( N1, N2 );

    common::Stencil2D<ValueType> stencil( 5 ); // take 5-point stencil

    // Stencil matrix is always square

    StencilStorage<ValueType> storage( grid, stencil );

    BOOST_CHECK_EQUAL( N1 * N2, storage.getNumRows() );
    BOOST_CHECK_EQUAL( N1 * N2, storage.getNumColumns() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( copyTest, ValueType, scai_numeric_test_types )
{
    const IndexType N1 = 5;
    const IndexType N2 = 8;
    common::Grid2D grid( N1, N2 );

    common::Stencil2D<ValueType> stencil( 5 ); // take 5-point stencil

    StencilStorage<ValueType> storage1( grid, stencil );
    StencilStorage<ValueType> storage2( storage1 );
    std::unique_ptr<StencilStorage<ValueType> > storage3( storage1.copy() );

    BOOST_CHECK_EQUAL( storage1.getNumRows(), storage2.getNumRows() );
    BOOST_CHECK_EQUAL( storage1.getNumRows(), storage3->getNumRows() );

    BOOST_CHECK_EQUAL( storage1.getNumColumns(), storage2.getNumColumns() );
    BOOST_CHECK_EQUAL( storage1.getNumColumns(), storage3->getNumColumns() );

    BOOST_CHECK_EQUAL( storage1.getStencil(), storage2.getStencil() );
    BOOST_CHECK_EQUAL( storage1.getStencil(), storage3->getStencil() );
    BOOST_CHECK_EQUAL( storage1.getGrid(), storage2.getGrid() );
    BOOST_CHECK_EQUAL( storage1.getGrid(), storage3->getGrid() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( convertTest, ValueType, scai_numeric_test_types )
{
    const IndexType N1 = 5;
    const IndexType N2 = 8;
    common::Grid2D grid( N1, N2 );

    common::Stencil2D<ValueType> stencil( 5 ); // take 5-point stencil

    StencilStorage<ValueType> stencilStorage( grid, stencil );

    auto csrStorage = convert<CSRStorage<ValueType>>( stencilStorage );

    // 5 entries for each grid element, but subract one for all border elements

    IndexType expectedNNZ = 5 * N1 * N2 - 2 * N1 - 2 * N2;
 
    BOOST_CHECK_EQUAL( expectedNNZ, csrStorage.getNumValues() );

    for ( IndexType i1 = 0; i1 < N1; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < N2; ++i2 )
        {
            IndexType pos = grid.linearPos( i1, i2 );

            // check diagonal element

            BOOST_CHECK_EQUAL( ValueType( 4 ), csrStorage.getValue( pos, pos ) );

            if ( i1 > 0 ) 
            {
                IndexType posUp = grid.linearPos( i1 - 1, i2 );
                BOOST_CHECK_EQUAL( ValueType( -1 ), csrStorage.getValue( pos, posUp ) );
            }

            if ( i1 < N1 - 1 ) 
            {
                IndexType posDown = grid.linearPos( i1 + 1, i2 );
                BOOST_CHECK_EQUAL( ValueType( -1 ), csrStorage.getValue( pos, posDown ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( transposeTest )
{
    typedef DefaultReal ValueType;

    for ( int b1 = 0; b1 < 2; b1++ )
    {
        for ( int b2 = 0; b2 < 2; b2++ )
        {
            const IndexType N1 = 6;
            const IndexType N2 = 8;

            common::Grid2D grid( N1, N2 );

            // transpose only if boundary conditions are same

            grid.setBorderType( 0, common::Grid::BorderType( b1 ), common::Grid::BorderType( b1 ) );
            grid.setBorderType( 1, common::Grid::BorderType( b2 ), common::Grid::BorderType( b2 ) );

            const int stencilData[9] = {  -1, -2, 1,
                                          -1, 4, 5,
                                           3, 2, 1   };

            common::Stencil2D<ValueType> stencil( 3, 3, stencilData );

            common::Stencil2D<ValueType> stencilT;
            stencilT.transpose( stencil );

            StencilStorage<ValueType> stencilStorage( grid, stencil );

            SCAI_LOG_INFO( logger, "Test stencil tranpose on this storage: " << stencilStorage )

            StencilStorage<ValueType> stencilStorageT( grid, stencilT );
        
            auto csrStorage = convert<CSRStorage<ValueType>>( stencilStorage );
            csrStorage.sortRows();
            csrStorage.compress();

            auto csrStorageT = convert<CSRStorage<ValueType>>( stencilStorageT );
            csrStorageT.sortRows();
            csrStorageT.compress();
        
            CSRStorage<ValueType> csrStorageT1;
            csrStorageT1.assignTranspose( csrStorage );

            BOOST_CHECK_EQUAL( csrStorageT.maxDiffNorm( csrStorageT1 ), 0 );
        }
    }
}
        
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getTest )
{
    typedef DefaultReal ValueType;

    const IndexType N1 = 6;
    const IndexType N2 = 8;

    common::Grid2D grid( N1, N2 );

    const int stencilData[9] = {  -1, -2, 1,
                                  -1, 4, 5,
                                   3, 2, 1   };

    common::Stencil2D<ValueType> stencil( 3, 3, stencilData );

    StencilStorage<ValueType> stencilStorage( grid, stencil );

    auto csrStorage = convert<CSRStorage<ValueType>>( stencilStorage );

    for ( IndexType i = 0; i < N1; i++ )
    {
        HArray<ValueType> row1;
        HArray<ValueType> row2;

        stencilStorage.getRow( row1, i );
        BOOST_CHECK_EQUAL( row1.size(), N1 * N2 );
        csrStorage.getRow( row2, i );

        BOOST_TEST( hostReadAccess( row1 ) == hostReadAccess( row2 ), per_element() );
    }

    for ( IndexType i = 0; i < N1; i++ )
    {
        for ( IndexType j = 0; j < N2; j++ )
        {
            BOOST_CHECK_EQUAL( stencilStorage.getValue( i, j ), csrStorage.getValue( i, j ) );
        }
    }
}
        
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/**
 * @file StencilUtilsTest.cpp
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
 * @brief Contains tests for the StencilKernel interface to be tested on different devices
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others

#include <scai/hmemo.hpp>
#include <scai/kregistry.hpp>
#include <scai/utilskernel.hpp>
#include <scai/sparsekernel/StencilUtils.hpp>
#include <scai/utilskernel/test/TestMacros.hpp>

#include <scai/hmemo/test/ContextFix.hpp>
#include <scai/common/Stencil.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using common::TypeTraits;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.StencilUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( setupTest )
{
    typedef DefaultReal ValueType;

    const IndexType n1 = 15;
    const IndexType n2 = 10;

    common::Grid2D grid( n1, n2 );

    common::Stencil2D<ValueType> stencil( 5 );  // 5 point-stencil

    hmemo::HArray<IndexType> gridInfo;
    hmemo::HArray<int> stencilInfo;
    hmemo::HArray<ValueType> stencilValues;

    StencilUtils::setup( gridInfo, stencilInfo, stencilValues, grid, stencil );

    hmemo::HArray<IndexType> expGridInfo( { n1, n2, n2, 1, 0, 0, 0, 0, 1, 1, 1, 1 } );
    int i_n2 = static_cast<IndexType>( n2 );  // narrowing is okay
    hmemo::HArray<int> expStencilInfo( { 0, -1, 0, 1, -1, 0, 1, 0, 0, 0, -1, 1, -i_n2, i_n2, 0 } );
    hmemo::HArray<ValueType> expStencilValues( { -1, -1, -1, -1, 4 } );

    SCAI_CHECK_EQUAL_ARRAY( gridInfo, expGridInfo );
    SCAI_CHECK_EQUAL_ARRAY( stencilInfo, expStencilInfo );
    SCAI_CHECK_EQUAL_ARRAY( stencilValues, expStencilValues );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( clearSetupTest )
{
    typedef DefaultReal ValueType;

    // zero-sized grid, identity stencil

    common::Grid1D grid( 0 );
    common::Stencil1D<ValueType> stencil( 1 );  // point-stencil

    hmemo::HArray<IndexType> gridInfo;
    hmemo::HArray<int> stencilInfo;
    hmemo::HArray<ValueType> stencilValues;

    StencilUtils::setup( gridInfo, stencilInfo, stencilValues, grid, stencil );

    hmemo::HArray<IndexType> expGridInfo( { 0, 1, 0, 0, 0, 0 } );
    hmemo::HArray<int> expStencilInfo( { 0, 0 } );
    hmemo::HArray<ValueType> expStencilValues( 1, ValueType( 0 ) );

    SCAI_CHECK_EQUAL_ARRAY( gridInfo, expGridInfo );
    SCAI_CHECK_EQUAL_ARRAY( stencilInfo, expStencilInfo );
    SCAI_CHECK_EQUAL_ARRAY( stencilValues, expStencilValues );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stencilGEMV3Test )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "stencilGEMV3Test<" << common::TypeTraits<ValueType>::id() << "> on " << *testContext )

    const IndexType n1 = 3;
    const IndexType n2 = 4;
    const IndexType n3 = 3;

    const common::Grid3D grid( n1, n2, n3 );

    const IndexType gridSize = grid.size();   // mult n1 * n2 * n3 only once

    common::Stencil3D<ValueType> stencil;

    stencil.addPoint(  0, 0, 0, 2 );
    stencil.addPoint(  1, 0, 0, -1 );   // right neighbor in 1st dim
    stencil.addPoint( -1, 0, 0, -1 );   // left neighbor in 1st dim

    const IndexType nPoints = stencil.nPoints();

    IndexType nDims = 3;

    hmemo::HArray<IndexType> gridInfo;
    hmemo::HArray<int> stencilInfo;
    hmemo::HArray<ValueType> stencilValues;

    StencilUtils::setup( gridInfo, stencilInfo, stencilValues, grid, stencil );

    HArray<IndexType> expGridInfo( { n1, n2, n3, n3 * n2, n3, 1, 
                                     0,  0,  0,  0,  0,  0,
                                     1,  1,  0,  0,  0,  0 } );

    SCAI_CHECK_EQUAL_ARRAY( gridInfo, expGridInfo );

    const HArray<ValueType> x( { 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1,
                                 1, 1, 2, 5, 3, 4, 1, 2, 3, 3, 0, 2,
                                 4, 1, 2, 3, 1, 0, 1, 2, 2, 3, 1, 1
                               } );

    BOOST_REQUIRE_EQUAL( x.size(), gridSize );

    HArray<ValueType> y1( gridSize, ValueType( 0 ) );   // result by hand
    HArray<ValueType> y2( gridSize, ValueType( 0 ) );   // result by kernel

    // apply stencil manually

    {
        WriteOnlyAccess<ValueType> wY( y2, gridSize );
        ReadAccess<ValueType> rX( x );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                for ( IndexType k = 0; k < n3; ++k )
                {
                    wY[ i * n2 * n3 + j * n3 + k ] = 2 * rX[ i * n2 * n3 + j * n3 + k ];

                    if ( i > 0 )
                    {
                        wY[ i * n2 * n3 + j * n3 + k ] -= rX[ ( i - 1 ) * n2 * n3 + j * n3 + k ];
                    }

                    if ( i < n1 - 1 )
                    {
                        wY[ i * n2 * n3 + j * n3 + k ] -= rX[ ( i + 1 ) * n2 * n3 + j * n3 + k ];
                    }
                }
            }
        }
    }

    // apply stencil by kernel

    const ValueType alpha = 1;
    const ValueType beta  = 1;

    StencilUtils::gemv( y1, alpha, x, beta, y1, 
                        gridSize, nDims, grid.sizes(), nPoints,
                        gridInfo, stencilInfo, stencilValues, false, testContext );

    BOOST_TEST( hostReadAccess( y1 ) == hostReadAccess( y2 ), per_element() );

    {
        ReadAccess<ValueType> rY1( y1 );
        ReadAccess<ValueType> rY2( y2 );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                for ( IndexType k = 0; k < n3; ++k )
                {
                    ValueType v1 = rY1[ i * n2 * n3 + j * n3 + k ];
                    ValueType v2 = rY2[ i * n2 * n3 + j * n3 + k ];

                    if ( v1 == v2 )
                    {
                        continue;
                    }

                    SCAI_LOG_ERROR( logger, "point (" << i << ", " << j << ", " << k << " ) : by kernel: " << v1 << ", by hand " << v2 )
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stencilGEMV4Test )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "stencilGEMV4Test<" << common::TypeTraits<ValueType>::id() << "> on " << *testContext )

    const IndexType n1 = 5;
    const IndexType n2 = 4;
    const IndexType n3 = 3;
    const IndexType n4 = 5;

    common::Grid4D grid( n1, n2, n3, n4 );
    grid.setBorderType( 3, common::BorderType::PERIODIC );

    const IndexType gridSize = grid.size();   // mult n1 * n2 * n3 only once

    common::Stencil4D<ValueType> stencil;

    stencil.addPoint(  0, 0, 0, 0, 3 );
    stencil.addPoint(  1, 0, 0, 0, -1 );   // right neighbor in 1st dim
    stencil.addPoint(  0, 0, 0, 1, 1 );   // right neighbor in last dim
    stencil.addPoint( -1, 0, 0, 0, -2 );   // left neighbor in 1st dim

    const IndexType nPoints = stencil.nPoints();

    IndexType nDims = 4;

    hmemo::HArray<IndexType> gridInfo;
    hmemo::HArray<int> stencilInfo;
    hmemo::HArray<ValueType> stencilValues;

    StencilUtils::setup( gridInfo, stencilInfo, stencilValues, grid, stencil );

    HArray<ValueType> x;
    HArrayUtils::setSequence<ValueType>( x, 1, 1, gridSize );

    BOOST_REQUIRE_EQUAL( x.size(), gridSize );

    HArray<ValueType> y1;   // result by kernel
    HArray<ValueType> y2;   // result by hand

    // apply stencil manually

    {   
        WriteOnlyAccess<ValueType> wY( y2, gridSize );
        ReadAccess<ValueType> rX( x );

        for ( IndexType i = 0; i < n1; ++i )
        {   
            for ( IndexType j = 0; j < n2; ++j )
            {   
                for ( IndexType k = 0; k < n3; ++k )
                {   
                    for ( IndexType m = 0; m < n4; ++m )
                    {   
                        IndexType pos = i * n2 * n3 * n4 + j * n3 * n4 + k * n4 + m;
                        
                        wY[ pos] = 3 * rX[ pos ];
                        
                        if ( i > 0 )
                        {   
                            wY[ pos ] -= 2 * rX[ pos - n2 * n3 * n4 ];
                        }
                        
                        if ( i < n1 - 1 )
                        {   
                            wY[ pos ] -= rX[ pos + n2 * n3 * n4 ];
                        }

                        if ( m < n4 - 1 )
                        {   
                            wY[ pos ] += rX[ pos + 1 ];
                        }
                        else
                        {
                            wY[ pos ] += rX[ pos + 1 - n4 ];  // Periodic
                        }
                    }
                }
            }
        }
    }

    // apply stencil by kernel

    const ValueType alpha = 1;
    const ValueType beta  = 0;

    StencilUtils::gemv( y1, alpha, x, beta, y1, 
                        gridSize, nDims, grid.sizes(), nPoints,
                        gridInfo, stencilInfo, stencilValues, false, testContext );

    BOOST_TEST( hostReadAccess( y1 ) == hostReadAccess( y2 ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

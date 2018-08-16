/**
 * @file StencilKernelTest.cpp
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
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using common::TypeTraits;
using common::Grid;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilKernelTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.StencilKernelTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilLocalTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilLocalSizes> stencilLocalSizes;
    LAMAKernel<StencilKernelTrait::stencilLocalCSR<ValueType> > stencilLocalCSR;

    ContextPtr loc = testContext;
    stencilLocalSizes.getSupportedContext( loc, stencilLocalCSR );

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType nDims = 1;
    IndexType gridSizes[] = { 100 };
    IndexType gridDistances[] = { 1 };
    Grid::BorderType gridBorders[] = { Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING };
    IndexType nPoints = 3;
    int stencilNodes[] = { -1, 0, 1 };
    // ValueType stencilValues[] = { -1, 2, -1 };

    {
        IndexType n = gridSizes[0];   // number of grid points
        WriteOnlyAccess<IndexType> wIA( csrIA, loc, n + 1 );
        stencilLocalSizes[loc]( wIA.get(), nDims, gridSizes, gridDistances, gridBorders, nPoints, stencilNodes );
        wIA.resize( n );
    }

    BOOST_CHECK_EQUAL( 3, csrIA[1] );   // inner point: all 3 stencil points are valid
    BOOST_CHECK_EQUAL( 2, csrIA[0] );   // left boundary, only 2 points
    BOOST_CHECK_EQUAL( 2, csrIA[99] );  // right boundary, only 2 points

    // now build the CSR offsets

    IndexType nnz = HArrayUtils::scan1( csrIA );

    {
        WriteOnlyAccess<IndexType> wJA( csrJA, loc, nnz );
        WriteOnlyAccess<ValueType> wValues( csrValues, loc, nnz );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilGEMV1Test, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilGEMV<ValueType> > stencilGEMV;

    ContextPtr loc = testContext;
    stencilGEMV.getSupportedContext( loc );

    const IndexType n1 = 8;

    IndexType nDims = 1;
    IndexType gridSizes[] = { n1 };
    IndexType gridDistances[] = { 1 };
    int stencilNodes[] = { 0, -1,  1, };
    ValueType stencilValues[] = { 2, -1, -1 };

    const ValueType xRaw[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    const IndexType n = sizeof( xRaw ) / sizeof( ValueType );
    const IndexType nPoints = sizeof( stencilValues ) / sizeof( ValueType );
    const IndexType nNodes = sizeof( stencilNodes ) / sizeof( int );

    BOOST_CHECK_EQUAL( n, n1 );
    BOOST_CHECK_EQUAL( nPoints * nDims, nNodes );
  
    // For one-dimensional stencil we check for all border types 

    for ( int b1 = 0; b1 < 2; b1++ )
    {
        for ( int b2 = 0; b2 < 2; b2++ )
        {
            Grid::BorderType gridBorders[] = { Grid::BorderType( b1 ), Grid::BorderType( b2 ) };

            HArray<ValueType> x ( n, xRaw );

            HArray<ValueType> y1 ( n, ValueType( 0 ) );   // result by kernel, must be initialized
            HArray<ValueType> y2;                         // result by manual evaluation

            // apply stencil manually

            {
                WriteOnlyAccess<ValueType> wY( y2, n );
                ReadAccess<ValueType> rX( x );

                for ( IndexType i = 0; i < n1; ++i )
                {
                    wY[ i ] = 2 * rX[ i ];

                    if ( i > 0 )
                    {
                        wY[ i ] -= rX[ i - 1  ];
                    }
                    else if ( gridBorders[0] == Grid::BORDER_PERIODIC )
                    {
                        wY[ i ] -= rX[ n1 - 1 ];
                    }

                    if ( i < n1 - 1 )
                    {
                        wY[ i ] -= rX[ i + 1 ];
                    }
                    else if ( gridBorders[1] == Grid::BORDER_PERIODIC )
                    {
                        wY[ i ] -= rX[ 0 ];
                    }
                }
            }

            // apply stencil by kernel

            IndexType width[] = { 1, 1 };
            int stencilOffset[] = { 0, -1, 1 };

            {
                SCAI_CONTEXT_ACCESS( loc );

                WriteAccess<ValueType> wY( y1, loc, n );
                ReadAccess<ValueType> rX( x, loc );

                ValueType alpha = 1;
                IndexType nPoints = 3;

                stencilGEMV[loc]( wY.get(), alpha, rX.get(), nDims, gridSizes, width, gridDistances, gridBorders,
                                  nPoints, stencilNodes, stencilValues, stencilOffset );
            }

            BOOST_TEST( hostReadAccess( y1 ) == hostReadAccess( y2 ), per_element() );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilGEMV2Test, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilGEMV<ValueType> > stencilGEMV;

    ContextPtr loc = testContext;
    stencilGEMV.getSupportedContext( loc );

    const IndexType n1 = 4;
    const IndexType n2 = 5;

    IndexType nDims = 2;
    IndexType gridSizes[] = { n1, n2 };
    IndexType gridDistances[] = { n2, 1 };

    // Each combination of border types should work correctly

    Grid::BorderType gridBorders[] = { Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING,
                                       Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING
                                     };
    int stencilNodes[] = { 0, 0, -1, 0, 1, 0, 0, -1, 0, 1 };
    ValueType stencilValues[] = { 6, -1, -1, -2, -2 };

    const ValueType xRaw[] = { 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1, 1, 1, 1, 2, 3, 1, 4, 1 };

    const IndexType n = sizeof( xRaw ) / sizeof( ValueType );
    const IndexType nPoints = sizeof( stencilValues ) / sizeof( ValueType );
    const IndexType nNodes = sizeof( stencilNodes ) / sizeof( int );

    BOOST_CHECK_EQUAL( n, gridSizes[0] * gridSizes[1] );
    BOOST_CHECK_EQUAL( nPoints * nDims, nNodes );

    HArray<ValueType> x ( n, xRaw );

    HArray<ValueType> y1 ( n, ValueType( 0 ) );   // result by kernel, must be initialized
    HArray<ValueType> y2;                         // result by hand

    // apply stencil manually

    {
        WriteOnlyAccess<ValueType> wY( y2, n );
        ReadAccess<ValueType> rX( x );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                wY[ i * n2 + j ] = 6 * rX[ i * n2 + j ];

                if ( i > 0 )
                {
                    wY[ i * n2 + j ] -= rX[ ( i - 1 ) * n2 + j ];
                }
                else if ( gridBorders[0] == Grid::BORDER_PERIODIC )
                {
                    wY[ i * n2 + j ] -= rX[ ( n1 - 1 ) * n2 + j ];
                }

                if ( i < n1 - 1 )
                {
                    wY[ i * n2 + j ] -= rX[ ( i + 1 ) * n2 + j ];
                }
                else if ( gridBorders[1] == Grid::BORDER_PERIODIC )
                {
                    wY[ i * n2 + j ] -= rX[ ( 0 ) * n2 + j ];
                }

                if ( j > 0 )
                {
                    wY[ i * n2 + j ] -= 2 * rX[ i * n2 + j - 1 ];
                }

                if ( j < n2 - 1 )
                {
                    wY[ i * n2 + j ] -= 2 * rX[ i * n2 + j + 1 ];
                }
            }
        }
    }

    // apply stencil by kernel

    IndexType width[] = { 1, 1, 1, 1 };

    int tmpOffset2 = n2;
    int tmpOffset1 = 1;
    int stencilOffset[] = { 0, -tmpOffset2, tmpOffset2, -tmpOffset1, tmpOffset1 };

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wY( y1, loc, n );
        ReadAccess<ValueType> rX( x, loc );

        ValueType alpha = 1;
        stencilGEMV[loc]( wY.get(), alpha, rX.get(), nDims, gridSizes, width, gridDistances, gridBorders,
                          nPoints, stencilNodes, stencilValues, stencilOffset );
    }

    BOOST_TEST( hostReadAccess( y1 ) == hostReadAccess( y2 ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilGEMV3Test, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilGEMV<ValueType> > stencilGEMV;

    ContextPtr loc = testContext;
    stencilGEMV.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "stencilGEMV3Test<" << common::TypeTraits<ValueType>::id() << "> on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 4;
    const IndexType n3 = 3;

    IndexType nDims = 3;
    IndexType gridSizes[] = { n1, n2, n3 };
    IndexType gridDistances[] = { n2 * n3, n3, 1 };
    Grid::BorderType gridBorders[] = { Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING,
                                       Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING,
                                       Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING
                                     };
    int stencilNodes[] = { 0, 0, 0, -1, 0, 0, 1, 0, 0 };
    ValueType stencilValues[] = { 2, -1, -1 };

    const ValueType xRaw[] = { 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1,
                               1, 1, 2, 5, 3, 4, 1, 2, 3, 3, 0, 2,
                               4, 1, 2, 3, 1, 0, 1, 2, 2, 3, 1, 1
                             };

    const IndexType n = sizeof( xRaw ) / sizeof( ValueType );
    const IndexType nPoints = sizeof( stencilValues ) / sizeof( ValueType );
    const IndexType nNodes = sizeof( stencilNodes ) / sizeof( int );

    BOOST_REQUIRE_EQUAL( n, gridSizes[0] * gridSizes[1] * gridSizes[2] );
    BOOST_REQUIRE_EQUAL( nPoints * nDims, nNodes );

    HArray<ValueType> x ( n, xRaw );

    HArray<ValueType> y1 ( n, ValueType( 0 ) );   // result by hand, must be initialized
    HArray<ValueType> y2;   // result by kernel

    // apply stencil manually

    {
        WriteOnlyAccess<ValueType> wY( y2, n );
        ReadAccess<ValueType> rX( x );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                for ( IndexType k = 0; k < n3; ++k )
                {
                    wY[ i * n2 * n3 + j * n3 + k ] = 2 * rX[ i * n2 * n3 + j * n3 + k ];

                    continue;

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

    IndexType width[] = { 1, 1, 0, 0, 1, 0 };

    int tmpOffset = n2 * n3;
    int stencilOffset[] = { 0, -tmpOffset, tmpOffset };

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wY( y1, loc, n );
        ReadAccess<ValueType> rX( x, loc );

        ValueType alpha = 1;
        const IndexType nPoints = 1;
        stencilGEMV[loc]( wY.get(), alpha, rX.get(), nDims, gridSizes, width, gridDistances, gridBorders,
                          nPoints, stencilNodes, stencilValues, stencilOffset );
    }

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

BOOST_AUTO_TEST_CASE_TEMPLATE( stencilGEMV4Test, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<StencilKernelTrait::stencilGEMV<ValueType> > stencilGEMV;

    ContextPtr loc = testContext;
    stencilGEMV.getSupportedContext( loc );

    const IndexType n1 = 3;
    const IndexType n2 = 2;
    const IndexType n3 = 3;
    const IndexType n4 = 2;

    IndexType nDims = 4;
    IndexType gridSizes[] = { n1, n2, n3, n4 };
    IndexType gridDistances[] = { n2* n3 * n4, n3 * n4, n4, 1 };
    Grid::BorderType gridBorders[] = { Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING,
                                       Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING,
                                       Grid::BORDER_ABSORBING, Grid::BORDER_ABSORBING
                                     };
    int stencilNodes[] = { 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0 };
    ValueType stencilValues[] = { 2, -1, -1 };

    const ValueType xRaw[] = { 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1,
                               1, 1, 2, 5, 3, 4, 1, 2, 3, 3, 0, 2,
                               4, 1, 2, 3, 1, 0, 1, 2, 2, 3, 1, 1
                             };

    const IndexType n = sizeof( xRaw ) / sizeof( ValueType );
    const IndexType nPoints = sizeof( stencilValues ) / sizeof( ValueType );
    const IndexType nNodes = sizeof( stencilNodes ) / sizeof( int );

    BOOST_REQUIRE_EQUAL( n, gridSizes[0] * gridSizes[1] * gridSizes[2] * gridSizes[3] );
    BOOST_REQUIRE_EQUAL( nPoints * nDims, nNodes );

    HArray<ValueType> x ( n, xRaw );

    HArray<ValueType> y1 ( n, ValueType( 0 ) );   // result by hand, must be initialized
    HArray<ValueType> y2;   // result by kernel

    // apply stencil manually

    {
        WriteOnlyAccess<ValueType> wY( y2, n );
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

                        wY[ pos] = 2 * rX[ pos ];

                        if ( i > 0 )
                        {
                            wY[ pos ] -= rX[ pos - n2 * n3 * n4 ];
                        }

                        if ( i < n1 - 1 )
                        {
                            wY[ pos ] -= rX[ pos + n2 * n3 * n4 ];
                        }
                    }
                }
            }
        }
    }

    // apply stencil by kernel

    IndexType width[] = { 1, 1, 0, 0, 0, 0, 0, 0 };

    int tmpOffset = n2 * n3 * n4;
    int stencilOffset[] = { 0, -tmpOffset, tmpOffset };

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wY( y1, loc, n );
        ReadAccess<ValueType> rX( x, loc );

        ValueType alpha = 1;
        stencilGEMV[loc]( wY.get(), alpha, rX.get(), nDims, gridSizes, width, gridDistances, gridBorders,
                          nPoints, stencilNodes, stencilValues, stencilOffset );
    }

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
                    for ( IndexType m = 0; m < n4; ++m )
                    {
                        IndexType pos = i * n2 * n3 * n4 + j * n3 * n4 + k * n4 + m;
                        ValueType v1 = rY1[ pos ];
                        ValueType v2 = rY2[ pos ];

                        if ( v1 == v2 )
                        {
                            continue;
                        }

                        SCAI_LOG_ERROR( logger, "point (" << i << ", " << j << ", " << k << ", " << m
                                        << " ) : by kernel: " << v1 << ", by hand " << v2 )
                    }
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

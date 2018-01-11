/**
 * @file StencilTest.cpp
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
 * @brief Test routines for Stencil classes.
 * @author Thomas Brandes
 * @date 06.05.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Stencil.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <memory>

using namespace scai;

using common::Stencil;
using common::Stencil1D;
using common::Stencil2D;
using common::Stencil3D;
using common::Stencil4D;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    for ( IndexType ndim = 1; ndim < 5; ndim++ )
    {
        Stencil<ValueType> stencil( ndim );
        BOOST_CHECK_EQUAL( ndim, stencil.nDims() );
        BOOST_CHECK_EQUAL( 0, stencil.nPoints() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setMatrixTest1, ValueType, scai_numeric_test_types )
{   
    // be careful, this test does not work for int stencils

    const double FD8[9] = { 0,
                            -5.0/7168.0, 
                             49.0/5120.0, 
                             245.0/3072.0, 
                             1225.0/1024.0, 
                            -1225.0/1024.0, 
                             245.0/3072.0,
                            -49.0/5120.0, 
                             5.0/7168.0 };

    Stencil1D<ValueType> stencil( 9, FD8 );

    // only 8 points are relevant 

    BOOST_CHECK_EQUAL( 8, stencil.nPoints() );

    IndexType width[2];

    stencil.getWidth( width );

    BOOST_CHECK_EQUAL( 3, width[0] );
    BOOST_CHECK_EQUAL( 4, width[1] );

    BOOST_REQUIRE_EQUAL( 8, stencil.getMatrixSize() );

    ValueType matrix[8];

    stencil.getMatrix( matrix );

    for ( IndexType i = 0; i < 8; ++i )
    {
        BOOST_CHECK_CLOSE( FD8[i+1], static_cast<double>( matrix[i] ), 0.1 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setMatrixTest2, ValueType, scai_array_test_types )
{   
    const int findEdges[25] = { 0,  0, -1,  0,  0,
                                0,  0, -1,  0,  0,
                                0,  0,  2,  0,  0,
                                0,  0,  0,  0,  0,
                                0,  0,  0,  0,  0,
                              };

    Stencil2D<ValueType> stencil( 5, 5, findEdges );

    // only 3 points are relevant 

    BOOST_CHECK_EQUAL( 3, stencil.nPoints() );

    IndexType width[4] = { invalidIndex, invalidIndex, invalidIndex, invalidIndex };

    stencil.getWidth( width );

    BOOST_CHECK_EQUAL( 2, width[0] );
    BOOST_CHECK_EQUAL( 0, width[1] );
    BOOST_CHECK_EQUAL( 0, width[2] );
    BOOST_CHECK_EQUAL( 0, width[3] );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setMatrixTest3, ValueType, scai_numeric_test_types )
{   
    const int stencilData[27] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, 26, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1
                                };

    Stencil3D<ValueType> stencil( 3, 3, 3, stencilData );

    // all points are relevant 

    BOOST_CHECK_EQUAL( 27, stencil.nPoints() );

    // and it is the same as the default 27 point stencil

    Stencil3D<ValueType> stencil27( 27 );

    BOOST_CHECK_EQUAL( stencil, stencil27 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stencil2Test )
{
    typedef float ValueType;

    Stencil1D<ValueType> stencil1( 3 );   // stencil in one-direction

    Stencil2D<ValueType> stencil2a( 5 );   
    Stencil2D<ValueType> stencil2b( stencil1, stencil1 );   

    BOOST_CHECK_EQUAL( stencil2a, stencil2b );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stencil3Test )
{
    typedef float ValueType;

    Stencil1D<ValueType> stencil1( 3 );   // stencil in one-direction

    Stencil3D<ValueType> stencil3a( 7 );
    Stencil3D<ValueType> stencil3b( stencil1, stencil1, stencil1 );

    BOOST_CHECK_EQUAL( stencil3a, stencil3b );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stencil4Test )
{
    typedef double ValueType;

    Stencil1D<ValueType> stencil1( 3 );   // stencil in one-direction

    Stencil4D<ValueType> stencil4a( 9 );   
    Stencil4D<ValueType> stencil4b( stencil1, stencil1, stencil1, stencil1 );

    BOOST_CHECK_EQUAL( stencil4a, stencil4b );

    Stencil4D<ValueType> stencil4c;

    stencil4c.addPoint( 0, 0, 0, 0, 8 );
    stencil4c.addPoint( 0, 0, 0, 1, -1 );
    stencil4c.addPoint( 0, 0, 0, -1, -1 );
    stencil4c.addPoint( 0, 0, 1, 0, -1 );
    stencil4c.addPoint( 0, 0, -1, 0,  -1 );
    stencil4c.addPoint( 0, 1, 0 , 0, -1 );
    stencil4c.addPoint( 0, -1, 0, 0,  -1 );

    // some points are misssing, so not equal 

    BOOST_CHECK( stencil4c != stencil4a );

    // add the missing two points

    stencil4c.addPoint( 1, 0, 0 , 0, -1 );
    stencil4c.addPoint( -1, 0, 0, 0,  -1 );
 
    // now we have the same stencil

    BOOST_CHECK( stencil4c == stencil4a );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validPointTest )
{
    typedef double ValueType;

    const IndexType nPoints = 7;

    Stencil3D<ValueType> stencil3( nPoints );

    const IndexType gridSizes[3] = { 4, 5, 4 };

    IndexType pos[3];

    std::unique_ptr<bool[]> valid( new bool[nPoints] );

    for ( IndexType i0 = 0; i0 < gridSizes[0]; ++i0 )
    {
        pos[0] = i0;

        for ( IndexType i1 = 0; i1 < gridSizes[1]; ++i1 )
        {
            pos[1] = i1;

            for ( IndexType i2 = 0; i2 < gridSizes[2]; ++i2 )
            {
                pos[2] = i2;
                
                const IndexType validPoints = stencil3.getValidPoints( valid.get(), gridSizes, pos );

                BOOST_CHECK( validPoints <= nPoints );

                IndexType cnt = 0;

                for ( IndexType k = 0; k < nPoints; ++k )
                {
                    if ( valid[k] ) 
                    { 
                        cnt++;
                    }
                    else
                    {
                        int p0, p1, p2;

                        ValueType v;
 
                        stencil3.getPoint( p0, p1, p2, v, k );

                        if ( false )
                        {
                            std::cout << "Point( " << i0 << ", " << i1 << ", " << i2 
                                      << " ) has no valid neighbor for stencil point " << k 
                                      << " = (" << p0 << ", " << p1 << ", " << p2 << " )" << std::endl;
                        }

                        if ( p0 == 1 )
                        {
                            BOOST_CHECK_EQUAL( i0 + 1, gridSizes[0] );
                        }
                        if ( p0 == -1 )
                        {
                            BOOST_CHECK_EQUAL( i0, IndexType( 0 ) );
                        }
                    }
                }

                BOOST_CHECK_EQUAL( cnt, validPoints );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeTest, ValueType, scai_array_test_types )
{   
    const int stencilData[9] = {  -1, -2, 1,
                                  -1, 4, 5,
                                   3, 2, 1   };

    // transposed stencil is given by point symmetry

    const int stencilDataT[9] = {   1,  2,  3,
                                    5,  4, -1,
                                    1, -2, -1  };

    Stencil2D<ValueType> stencil( 3, 3, stencilData );
    Stencil2D<ValueType> stencilT1( 3, 3, stencilDataT );
    Stencil2D<ValueType> stencilT2;
    stencilT2.transpose( stencil );

    BOOST_CHECK_EQUAL( stencilT1, stencilT2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_array_test_types )
{
    const int stencilData[9] = {  -1, -2, 1,
                                  -1, 4, 3,
                                   3, 2, 1   };

    // scaled stencil is given by scaling each value

    const int stencilDataS[9] = {  -2, -4, 2,
                                   -2, 8, 6,
                                    6, 4, 2   };

    Stencil2D<ValueType> stencilS1( 3, 3, stencilDataS );
    Stencil2D<ValueType> stencilS2( 3, 3, stencilData );
    stencilS2.scale( 2 );

    BOOST_CHECK_EQUAL( stencilS1, stencilS2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeTest )
{
    std::ostringstream f;

    Stencil3D<int> stencil( 27 );

    f << stencil;

    const std::string& fstr = f.str();

    BOOST_CHECK( fstr.length() > 7 );

    BOOST_CHECK( fstr.find( "27" ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();


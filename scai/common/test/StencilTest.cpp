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

using scai::common::Stencil;
using scai::common::Stencil1D;
using scai::common::Stencil2D;
using scai::common::Stencil3D;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_array_test_types )
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

    IndexType lb[1];
    IndexType ub[1];

    stencil.getWidth( lb, ub );

    BOOST_CHECK_EQUAL( 3, lb[0] );
    BOOST_CHECK_EQUAL( 4, ub[0] );

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

    IndexType lb[2];
    IndexType ub[2];

    stencil.getWidth( lb, ub );

    BOOST_CHECK_EQUAL( 2, lb[0] );
    BOOST_CHECK_EQUAL( 0, lb[1] );
    BOOST_CHECK_EQUAL( 0, ub[0] );
    BOOST_CHECK_EQUAL( 0, ub[1] );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setMatrixTest3, ValueType, scai_array_test_types )
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


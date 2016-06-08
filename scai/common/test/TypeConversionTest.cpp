/**
 * @file test/TypeConversionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Test for type conversions
 * @author Eric Schricker
 * @date 11.04.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace common;

#ifdef SCAI_COMPLEX_SUPPORTED

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TypeConversionTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Complex2ScalarTest )
{
    ComplexFloat c;
    ComplexDouble z;
    ComplexLongDouble lz;
    float f;
    double d;
    long double l;
    /*
     * Complex Float
     */
    c = ComplexFloat( 2.3, 3 );
    f = static_cast<float>( c );
    BOOST_CHECK_CLOSE( f, 2.3, 1e-5 );
    d = static_cast<double>( c );
    BOOST_CHECK_CLOSE( d, 2.3, 1e-5 );
    l = static_cast<long double>( c );
    BOOST_CHECK_CLOSE( l, 2.3, 1e-5 );
    /*
     * ComplexDouble
     */
    z = ComplexDouble( 2.3, 2.0f );
    f = static_cast<float>( z );
    BOOST_CHECK_CLOSE( f, 2.3, 1e-5 );
    d = static_cast<double>( z );
    BOOST_CHECK_CLOSE( d, 2.3, 1e-312 );
    l = static_cast<long double>( z );
    BOOST_CHECK_CLOSE( l, 2.3, 1e-312 );
    /*
     * ComplexLongDouble
     */
    lz = ComplexLongDouble( 2.3, 2.0f );
    f = static_cast<float>( lz );
    BOOST_CHECK_CLOSE( f, 2.3, 1e-5 );
    d = static_cast<double>( lz );
    BOOST_CHECK_CLOSE( d, 2.3, 1e-312 );
    l = static_cast<long double>( lz );
    BOOST_CHECK_CLOSE( l, 2.3, 1e-320 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Scalar2ComplexTest )
{
    float f;
    double d;
    long double l;
    ComplexFloat c;
    ComplexDouble z;
    ComplexLongDouble lz;
    /*
     * float
     */
    f = 3.31f;
    c = static_cast<ComplexFloat>( f );
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    c = f;
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    z = static_cast<ComplexDouble>( f );
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-5 );
    z = f;
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-5 );
    lz = static_cast<ComplexDouble>( f );
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-5 );
    lz = f;
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-5 );
    /*
     * double
     */
    d = 3.31;
    c = static_cast<ComplexFloat>( d );
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    c = d;
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    z = static_cast<ComplexDouble>( d );
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-312 );
    z = d;
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-312 );
    lz = static_cast<ComplexDouble>( d );
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-312 );
    lz = d;
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-312 );
    /*
     * long double
     */
    l = 3.31;
    c = static_cast<ComplexFloat>( l );
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    c = l;
    BOOST_CHECK_CLOSE( c.real(), 3.31, 1e-5 );
    BOOST_CHECK_CLOSE( c.imag(), 0.0, 1e-5 );
    z = static_cast<ComplexDouble>( l );
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-312 );
    z = l;
    BOOST_CHECK_CLOSE( z.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( z.imag(), 0.0, 1e-312 );
    lz = static_cast<ComplexDouble>( l );
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-312 );
    lz = l;
    BOOST_CHECK_CLOSE( lz.real(), 3.31, 1e-312 );
    BOOST_CHECK_CLOSE( lz.imag(), 0.0, 1e-312 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Complex2ComplexTest )
{
    ComplexFloat c_i, c_o;
    ComplexDouble z_i, z_o;
    ComplexLongDouble lz_i, lz_o;
    /*
     * ComplexFloat
     */
    c_i = ComplexFloat( -3, 2.0 );
    c_o = static_cast<ComplexFloat>( c_i );
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    c_o = c_i;
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    z_o = static_cast<ComplexDouble>( c_i );
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-5 );
    z_o = c_i;
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-5 );
    lz_o = static_cast<ComplexLongDouble>( c_i );
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-5 );
    lz_o = c_i;
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-5 );
    /*
     * ComplexDouble
     */
    z_i = ComplexDouble( -3, 2.0 );
    c_o = static_cast<ComplexFloat>( z_i );
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    c_o = z_i;
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    z_o = static_cast<ComplexDouble>( z_i );
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-312 );
    z_o = z_i;
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-312 );
    lz_o = static_cast<ComplexLongDouble>( z_i );
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-312 );
    lz_o = z_i;
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-312 );
    /*
     * ComplexLongDouble
     */
    lz_i = ComplexLongDouble( -3, 2.0 );
    c_o = static_cast<ComplexFloat>( lz_i );
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    c_o = lz_i;
    BOOST_CHECK_CLOSE( c_o.real(), -3, 1e-5 );
    BOOST_CHECK_CLOSE( c_o.imag(), 2, 1e-5 );
    z_o = static_cast<ComplexDouble>( lz_i );
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-312 );
    z_o = lz_i;
    BOOST_CHECK_CLOSE( z_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( z_o.imag(), 2, 1e-312 );
    lz_o = static_cast<ComplexLongDouble>( lz_i );
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-312 );
    lz_o = lz_i;
    BOOST_CHECK_CLOSE( lz_o.real(), -3, 1e-312 );
    BOOST_CHECK_CLOSE( lz_o.imag(), 2, 1e-312 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

#endif


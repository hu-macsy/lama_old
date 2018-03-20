/**
 * @file FFTTest.cpp
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
 * @brief Tests for FFT operations on arrays
 * @author Eric Schricker
 * @date 22.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/FFTKernelTrait.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Constants.hpp>
#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ContextAccess.hpp>

using namespace scai;
using namespace utilskernel;
using namespace common;
using namespace hmemo;

extern ContextPtr testContext;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( FFTTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.FFTTest" )

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<double> scai_fft_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( paddedForward1D, ValueType, scai_fft_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    static LAMAKernel<FFTKernelTrait::paddedForward1D<ValueType> > fft;

    fft.getSupportedContext( loc );

    const IndexType n = 2;
    const IndexType npad = 8;

    // Init Arrays

    HArray<ValueType> in( { 0.2, 0.16 } );

    HArray<Complex<ValueType> > out;

    // FFT Call
    {
        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<ValueType> rIn( in, loc );
        WriteOnlyAccess<Complex<ValueType> > rOut( out, loc, npad );
        fft[loc](n, npad, rIn.get(), rOut.get());
    }

    ValueType eps = 0.00001;

    for ( IndexType i = 0; i < npad; ++i )
    {
        Complex<ValueType> v = out[i];
        SCAI_LOG_ERROR( logger, "rout[" << i << "] = " << v )
    }

    // Test
    {
        ReadAccess<Complex<ValueType> > rOut( out );

        BOOST_CHECK( common::Math::abs( rOut[0] - Complex<ValueType>( 0.2 + 0.16 ) ) < eps );
        BOOST_CHECK( common::Math::abs( rOut[npad-2] - Complex<ValueType>( 0.2, 0.16 ) ) < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( paddedBackward1D, ValueType, scai_fft_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    static LAMAKernel<FFTKernelTrait::paddedBackward1D<ValueType> > fft;

    fft.getSupportedContext( loc );

    const IndexType n = 2;
    const IndexType npad = 8;

    // Init Arrays

    HArray<ValueType> in( { 0.2, 0.16 } );

    HArray<Complex<ValueType> > out;

    // FFT Call
    {
        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<ValueType> rIn( in, loc );
        WriteOnlyAccess<Complex<ValueType> > rOut( out, loc, npad );
        fft[loc](n, npad, rIn.get(), rOut.get());
    }

    ValueType eps = 0.00001;

    for ( IndexType i = 0; i < npad; ++i )
    {
        Complex<ValueType> v = out[i];
        SCAI_LOG_ERROR( logger, "rout[" << i << "] = " << v )
    }

    // Test
    {
        ReadAccess<Complex<ValueType> > rOut( out );

        BOOST_CHECK( common::Math::abs( rOut[0] - Complex<ValueType>( 0.2 + 0.16 ) ) < eps );
        BOOST_CHECK( common::Math::abs( rOut[npad-2] - Complex<ValueType>( 0.2, -0.16 ) ) < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

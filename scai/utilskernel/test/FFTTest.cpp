/**
 * @file FFTTest.cpp
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
 * @brief Tests for FFT operations on arrays
 * @author Eric Schricker
 * @date 22.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/FFTKernelTrait.hpp>
#include <scai/utilskernel/FFTUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Constants.hpp>

#include <scai/hmemo.hpp>

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

#ifdef SCAI_COMPLEX_SUPPORTED

// only use those ValueType where also ComplexType<RealType<ValueType>> is supported

typedef boost::mpl::list<SCAI_FFT_TYPES_HOST> scai_fft_test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftForwardTest, ValueType, scai_fft_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    typedef common::Complex<RealType<ValueType>> ComplexType;

    // Init Arrays

    HArray<ComplexType> data( { 0.2, 0.16, 0, 0, 0, 0, 0, 0 } );

    IndexType n = 8;

    FFTUtils::fftcall<ValueType>( data, 1, n, 3, 1, loc );

    RealType<ValueType> eps = 0.00001;

    {
        ReadAccess<ComplexType> rOut( data );

        BOOST_CHECK( common::Math::abs( rOut[0] - ComplexType( 0.2 + 0.16 ) ) < eps );
        BOOST_CHECK( common::Math::abs( rOut[n-2] - ComplexType( 0.2, 0.16 ) ) < eps );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftBackwardTest, ValueType, scai_fft_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    // Init Arrays

    typedef common::Complex<RealType<ValueType>> ComplexType;

    HArray<ComplexType> data( { 0.2, 0.16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  } );

    SCAI_LOG_DEBUG( logger, "in = " << common::Wrapper<hmemo::HostReadAccess<ComplexType>>( hmemo::hostReadAccess( data ) ) )

    IndexType n = data.size();

    FFTUtils::fftcall<ValueType>( data, 1, n, 3, -1, loc );

    SCAI_LOG_DEBUG( logger, "out = " << common::Wrapper<hmemo::HostReadAccess<ComplexType>>( hmemo::hostReadAccess( data ) ) )

    RealType<ValueType> eps = 0.00001;

    {
        ReadAccess<ComplexType> rOut( data );

        BOOST_CHECK( common::Math::abs( rOut[0] - ComplexType( 0.2 + 0.16 ) ) < eps );
        BOOST_CHECK( common::Math::abs( rOut[n-2] - ComplexType( 0.2, -0.16 ) ) < eps );
    }
}

/* --------------------------------------------------------------------- */

/* Discrete fourier transform ( inefficient straight-forward implementation ) */

template<typename ValueType>
static void discreteFourierTransform( Complex<ValueType> x2[], const Complex<ValueType> x1[], IndexType n, int dir )
{   
    ValueType PI_2 = static_cast<ValueType>( 2 * 3.14159265358979 );

    for ( IndexType i = 0; i < n; i++ )
    {
        x2[i] = 0;

        ValueType arg = - dir * PI_2 * ValueType( i ) / ValueType( n );

        for ( IndexType k = 0; k < n; k++ )
        {
            ValueType cosarg = std::cos( k * arg );
            ValueType sinarg = std::sin( k * arg );
            x2[i] += x1[k] * Complex<ValueType>( cosarg , sinarg );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTest, ValueType, scai_fft_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    typedef common::Complex<RealType<ValueType>> ComplexType;

    for ( int dir = -1; dir <=1; dir += 2 )
    {
        const IndexType m = 6;
        const IndexType n = 1 << m; 

        HArray<ComplexType> in( n );

        // generate random numbers between -1 and 1

        HArrayUtils::fillRandom( in, 2, 1.0f, loc );

        HArrayUtils::compute( in, in, common::BinaryOp::SUB, ComplexType( 1, 1 ) );

        HArray<ComplexType> outFFT( in );

        IndexType k = 1;

        FFTUtils::fftcall<ValueType>( outFFT, k, n, m, dir, loc );

        BOOST_REQUIRE_EQUAL( outFFT.size(), n );

        HArray<ComplexType> outDFT;

        {
            auto rIn = hostReadAccess( in );

            for ( IndexType i = 0; i < n; ++i )
            {
                // std::cout << "in[" << i << "] = " << rIn[i] << std::endl;
            }

            auto wOut = hostWriteOnlyAccess( outDFT, n );
            discreteFourierTransform( wOut.get(), rIn.get(), n, dir );
        }

        {
            auto rOut1 = hostReadAccess( outFFT );
            auto rOut2 = hostReadAccess( outDFT );

            RealType<ValueType> eps = common::TypeTraits<ComplexType>::small() * m;

            for ( IndexType i = 0; i < n; ++i )
            {
                auto diff = common::Math::abs( rOut1[i] - rOut2[i] );
    
                // std::cout << "FFT[ " << i << " ] = " << rOut1[i] << ", DFT[ " << i << " ] = " << rOut2[i] 
                //          << ", diff = " << diff << ", eps = " << eps << std::endl;

                BOOST_CHECK( diff < eps );
            }
        }
    }
}

#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

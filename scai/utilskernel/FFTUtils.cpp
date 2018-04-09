/**
 * @file FFTUtils.cpp
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
 * @brief Fast fourier transformation for heterorgeneous arrays
 * @author Thomas Brandes
 * @date 21.03.2018
 */

#include <scai/utilskernel/FFTUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/FFTKernelTrait.hpp>
#include <scai/utilskernel/SectionKernelTrait.hpp>

#include <scai/common/macros/loop.hpp>

namespace scai
{

using namespace hmemo;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( FFTUtils::logger, "FFTUtils" )

#ifdef SCAI_COMPLEX_SUPPORTED

using common::Complex;

/* --------------------------------------------------------------------------- */

static void pow2( IndexType& m, IndexType& n2, const IndexType n )
{
    m  = 0;
    n2 = 1;

    while ( n2 < n )
    {
        m  += 1;
        n2 = n2 << 1;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTUtils::fft( 
    HArray<Complex<RealType<ValueType>>>& result, 
    const HArray<ValueType>& x, 
    const IndexType n,
    int direction,
    const ContextPtr context )
{
    typedef Complex<RealType<ValueType>> FFTType;

    static LAMAKernel<UtilKernelTrait::set<FFTType, ValueType>> set;
    static LAMAKernel<UtilKernelTrait::setVal<FFTType>> setVal;
    static LAMAKernel<FFTKernelTrait::fft<RealType<ValueType>>> fft;

    ContextPtr loc = context;

    if ( !loc )
    {
        loc = x.getValidContext();
    }

    fft.getSupportedContext( loc, set, setVal );

    IndexType n2;
    IndexType m;
 
    pow2( m, n2, n );

    SCAI_LOG_INFO( logger, "fft<" << common::TypeTraits<ValueType>::id() 
                    << ">( array[" << n << " -> " << n2 << " = 2 ** " << m << "], dir = " << direction
                    << " @ ctx = " << *loc ) 

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<FFTType> wResult( result, loc, n2 );

    IndexType nx = std::min( n, x.size() );

    SCAI_LOG_DEBUG( logger, "copy x[" << nx << "] and pad to " << n2 )
    set[loc]( wResult.get(), rX.get(), x.size(), common::BinaryOp::COPY );
    setVal[loc]( wResult.get() + nx, n2 - nx, FFTType( 0 ), common::BinaryOp::COPY );

    fft[loc]( wResult.get(), 1, n2, m, direction );

    wResult.resize( n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTUtils::fft_many( 
    HArray<Complex<RealType<ValueType>>>& result, 
    const HArray<ValueType>& x, 
    const IndexType many,
    const IndexType ncols,
    int direction,
    const ContextPtr context )
{
    typedef Complex<RealType<ValueType>> FFTType;

    IndexType ncols2;
    IndexType m;
 
    pow2( m, ncols2, ncols );

    SCAI_ASSERT_EQ_ERROR( ncols2, ncols, "not supported yet" )

    SCAI_LOG_INFO( logger, "fft_many<" << common::TypeTraits<ValueType>::id() 
                     << ">( array[ " << many << " x " << ncols << "], " << ncols2 << " = 2 ** " << m << ", dir = " << direction )
 
    // copy input array x into result array and pad the rows with 0

    if ( static_cast<const _HArray*>( &result ) == static_cast<const _HArray*>( &x ) )
    {
        // we have an alias on input and output array

        SCAI_ASSERT_EQ_ERROR( x.size(), many * ncols2, "Size of aliased input/output not legal" )
    }
    else
    {
        IndexType nColsX = x.size() / many;
        SCAI_ASSERT_EQ_ERROR( nColsX * many, x.size(), "size of input array is not multiple of #many = " << many )

        HArrayUtils::setSameValue( result, many * ncols2, FFTType( 0 ), context );

        IndexType sizes[2] = { many, nColsX };
        IndexType sdist[2] = { nColsX, 1 };
        IndexType tdist[2] = { ncols2, 1 };

        static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::unaryOp<FFTType, ValueType> > unaryOp;

        ContextPtr loc = context ? context : x.getValidContext();

        unaryOp.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<ValueType> rX( x, loc );
        WriteAccess<FFTType> wResult( result, loc );

        unaryOp[loc]( wResult.get(), 2, sizes, tdist, rX.get(), sdist, common::UnaryOp::COPY );
    }

    // for ( IndexType i = 0; i < result.size(); ++i )
    // {
    //     FFTType v = result[i];
    //     std::cout << "before fft_n: result[" << i << "] = " << v << std::endl;
    // }

    {
        static LAMAKernel<FFTKernelTrait::fft<RealType<ValueType>>> fft;

        ContextPtr loc = context ? context : x.getValidContext();

        fft.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<FFTType> wResult( result, loc );
        fft[loc]( wResult.get(), many, ncols2, m, direction );
    }

    // for ( IndexType i = 0; i < result.size(); ++i )
    // {
    //    FFTType v = result[i];
    //    std::cout << "after fft_n: result[" << i << "] = " << v << std::endl;
    // }

}

/* --------------------------------------------------------------------------- */

#define FFTUTILS_SPECIFIER( ValueType )                     \
    template void FFTUtils::fft<ValueType>(                 \
        hmemo::HArray<Complex<RealType<ValueType>>>&,       \
        const hmemo::HArray<ValueType>&,                    \
        const IndexType n,                                  \
        const int direction,                                \
        hmemo::ContextPtr);                                 \
    template void FFTUtils::fft_many<ValueType>(            \
        hmemo::HArray<Complex<RealType<ValueType>>>&,       \
        const hmemo::HArray<ValueType>&,                    \
        const IndexType,                                    \
        const IndexType,                                    \
        const int,                                          \
        hmemo::ContextPtr);                    


    SCAI_COMMON_LOOP( FFTUTILS_SPECIFIER, SCAI_FFT_TYPES_HOST )

#undef FFTUTILS_SPECIFIER

#endif

} /* end namespace utilskernel */

} /* end namespace scai */


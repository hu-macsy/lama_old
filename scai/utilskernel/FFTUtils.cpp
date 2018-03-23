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

#include <scai/common/macros/loop.hpp>

namespace scai
{

using common::Complex;
using namespace hmemo;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( FFTUtils::logger, "FFTUtils" )

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

    SCAI_LOG_ERROR( logger, "fft<" << common::TypeTraits<ValueType>::id() 
                     << ">( array[" << n << "], " << n2 << " = 2 ** " << m << ", dir = " << direction )

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rX( x );
    WriteOnlyAccess<FFTType> wResult( result, loc, n2 );

    IndexType nx = std::min( n, x.size() );

    SCAI_LOG_ERROR( logger, "copy x[" << nx << "] and pad to " << n2 )
    set[loc]( wResult.get(), rX.get(), x.size(), common::BinaryOp::COPY );
    setVal[loc]( wResult.get() + nx, n2 - nx, FFTType( 0 ), common::BinaryOp::COPY );

    for ( IndexType i = 0; i < n2; ++i )
    {
        std::cout << "result[" << i << "] = " << wResult[i] << std::endl;
    }

    fft[loc]( wResult.get(), n2, m, direction );

    for ( IndexType i = 0; i < n2; ++i )
    {
        std::cout << "result[" << i << "] = " << wResult[i] << std::endl;
    }

    wResult.resize( n );
}

#define FFTUTILS_SPECIFIER( ValueType )                     \
    template void FFTUtils::fft<ValueType>(                 \
        hmemo::HArray<Complex<RealType<ValueType>>>&,       \
        const hmemo::HArray<ValueType>&,                    \
        const IndexType n,                                  \
        const int direction,                                \
        hmemo::ContextPtr);                    

SCAI_COMMON_LOOP( FFTUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef FFTUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */


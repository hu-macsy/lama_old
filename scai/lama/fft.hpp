/**
 * @file fft.hpp
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
 * @brief Fast fourier transformation for vector and matrices
 * @author Thomas Brandes
 * @date 21.03.2018
 */
#pragma once

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/utilskernel/FFTUtils.hpp>

namespace scai
{

namespace lama
{

/** Compute the discrete fourier transform of a dense vector using a fast fourier transform
 *
 *  @param[out] result is the result vector
 *  @param[in]  x  is the input vector
 *  @param[in]  n  padding length (n > x.size()) or truncate length , optional
 *
 */

template<typename ValueType>
void fft( 
    DenseVector<common::Complex<RealType<ValueType>>>& result, 
    const DenseVector<ValueType>& x, 
    const IndexType n = invalidIndex )
{
    // result array for fft is alway complex

    typedef common::Complex<RealType<ValueType>> FFTType;

    SCAI_ASSERT_ERROR( x.getDistribution().isReplicated(), "fft only on replicated vectors" )

    hmemo::ContextPtr ctx = result.getContextPtr();   // preferred location for FFT

    dmemo::DistributionPtr dummyDist;
    hmemo::HArray<FFTType> resultArray;

    result.splitUp( dummyDist, resultArray );    // reuse allocated array of result vector

    const IndexType size = n == invalidIndex ? x.size() : n;

    int direction = 1;

    utilskernel::FFTUtils::fft<ValueType>( resultArray, x.getLocalValues(), size, direction, result.getContextPtr() );

    result = DenseVector<FFTType>( std::move( resultArray ), ctx );
}

/**
template<typename ValueType>
void fft( 
    DenseMatrix<ComplexType<ValueType>>& result, 
    const DenseMatrix<ValueType>& x, 
    const IndexType dim,
    const IndexType n = invalidIndex )
{
    if ( dim == 0 )
    {
        // fft along the rows is done in parallel
        IndexType nLocalRows = x.getRowDistribution().getLocalSize();
        utilskernel::fft( nLocalRows, result.getLocalValues(), x.getLocalValues(), result.getContextPtr() );
    }
    else if ( dim == 1 )
    {
        // fft along the columns we do on the transpose
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal dim argument " << dim << ", must be 0 or 1" )
    }
}

template<typename ValueType>
void fft2( 
    DenseMatrix<ComplexType<ValueType>>& result, 
    const DenseMatrix<ValueType>& x, 
    const IndexType m = x.getNumRows(),
    const IndexType n = x.getNumColumns() )
{
    fft( result, x, 0, m );
    fft( result, result, 1, n );
}

**/


} /* end namespace lama */

} /* end namespace scai */


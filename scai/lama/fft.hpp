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

#ifdef SCAI_COMPLEX_SUPPORTED

/** Compute the discrete fourier transform of a dense vector using a fast fourier transform
 *
 *  @param[out] result is the result vector, length will be n 
 *  @param[in]  x  is the input vector
 *  @param[in]  n  padding length (n > x.size()) or truncate length , optional
 *
 *  - This operation is not available for distributed vectors
 *  - The result vector is always a complex vector while the input vector can be real or complex.
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

    utilskernel::FFTUtils::fft1D<ValueType>( resultArray, x.getLocalValues(), size, ctx );

    result = DenseVector<FFTType>( std::move( resultArray ), ctx );
}

/** Apply the inverse fast Fourier transform  to a dense vector 
 *
 *  @param[out] result is the result vector, length will be n 
 *  @param[in]  x  is the input vector, must 
 *  @param[in]  n  padding length (n > x.size()) or truncate length , optional
 *
 *  - This operation is not available for distributed vectors
 *  - The result vector is always a complex vector while the input vector can be real or complex.
 */
template<typename ValueType>
void ifft(
    DenseVector<ValueType>& result,
    const DenseVector<common::Complex<RealType<ValueType>>>& x,
    const IndexType n = invalidIndex )
{
    // result array for fft is alway complex

    SCAI_ASSERT_ERROR( x.getDistribution().isReplicated(), "fft only on replicated vectors" )

    hmemo::ContextPtr ctx = result.getContextPtr();   // preferred location for FFT

    dmemo::DistributionPtr dummyDist;
    hmemo::HArray<ValueType> resultArray;

    result.splitUp( dummyDist, resultArray );    // reuse allocated array of result vector

    const IndexType size = n == invalidIndex ? x.size() : n;

    utilskernel::FFTUtils::ifft1D<ValueType>( resultArray, x.getLocalValues(), size, result.getContextPtr() );

    result = DenseVector<ValueType>( std::move( resultArray ), ctx );
}

/* --------------------------------------------------------------------- */

/** 
 *  @brief Apply FFT on a two-dimensional matrix (in-place)
 *
 *  - For dim = 0 the FFT is applied to each column
 *  - For dim = 1 the FFT is applied to each row 
 *  - The matrix size in the dimension where the FFT is applied must be a power of 2
 *  - FFT can only be applied to complex matrices
 *  - use resize function to fill up, convert or truncate matrices
 */

template<typename ValueType>
void fftRows( 
    DenseMatrix<common::Complex<RealType<ValueType>>>& data,
    int direction )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    IndexType numColumns = data.getNumColumns();  // size of each vector to which FFT is applied

    pow2( m, n2, numColumns );
 
    SCAI_ASSERT_EQ_ERROR( n2, numColumns, "not power of 2" )

    if ( !data.getColDistribution().isReplicated() )
    {
        IndexType numColumns = data.getNumColumns();

        DistributionPtr repColumns = std::make_shared<NoDistribution( data.getNumColumns() );
        DistributionPtr distColumns = data.getColDistribuitonPtr();
        data.redistribute( data.getRowDistibutionPtr(), repColumns );
        fft( data, dim );
        data.redistribute( data.getRowDistibutionPtr(), distColumns );
    }
    else
    {
        // fft  for each row

        IndexType nLocalRows = data.getRowDistribution().getLocalSize();

        utilskernel::FFTUtils::fftcall( localData, nLocalRows, numColumns, m, direction, data.getContextPtr() );
    }
}

template<typename ValueType>
void fftCols(
    DenseMatrix<common::Complex<RealType<ValueType>>>& data,
    int direction )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    // transpose the matrix, apply FFT row-wise and transpose back

    if ( !data.getRowDistribution().isReplicated() && data.getColDistribution().isReplicated() )
    {
        dmemo::DistributionPtr saveColDistributionPtr = data.getColDistributionPtr();
        dmemo::DistributionPtr blockDist = std::make_shared<dmemo::BlockDistribution>( numColumns );
        data.redistribute( data.getRowDistributionPtr(), blockDist );
        fftCols( data, direction );
        data.redistribute( data.getRowDistribuitonPtr(), savedColDistributionPtr );
    }
    else if ( data.getRowDistribution().isReplicated() && !data.getColDistribution().isReplicated() )
    {
        dmemo::DistributionPtr saveRowDistributionPtr = data.getRowDistributionPtr();
        dmemo::DistributionPtr blockDist = std::make_shared<dmemo::BlockDistribution>( numRows );
        data.redistribute( blockDist, data.getColDistributionPtr() );
        fftCols( data, direction );
        data.redistribute( savedRowDistributionPtr, data.getColDistributionPtr() );
    }
    else
    {   // now transpose works fine, we do it in place

        data.assignTranspose( data );
        fftRows( data, direction );
        data.assignTranspose( data );
    }
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
void fft( 
    DenseMatrix<common::Complex<RealType<ValueType>>>& data,
    const IndexType dim )
{
   if ( dim == 0 )
   {
       fftCols( data, 1 ); // treat each column as a vector
   }
   else if ( dim == 1 )
   {
       fftRows( data, 1 ); // treat each row as a vector
   }
}

template<typename ValueType>
void ifft( 
    DenseMatrix<common::Complex<RealType<ValueType>>>& data,
    const IndexType dim )
{
   if ( dim == 0 )
   {
       fftCols( data, -1 ); // treat each column as a vector
   }
   else if ( dim == 1 )
   {
       fftRows( data, -1 ); // treat each row as a vector
   }
}

template<typename ValueType>
void fft(
    DenseMatrix<common::Complex<RealType<ValueType>>>& data )
{
    fftRows( data, -1 );   // fast backward FFT alon rows
    fftCols( data, -1 );   // fast backward FFT along columns
}

template<typename ValueType>
void ifft(
    DenseMatrix<common::Complex<RealType<ValueType>>>& data )
{
    fftRows( data, -1 );   // fast backward FFT alon rows
    fftCols( data, -1 );   // fast backward FFT along columns
}

#endif

} /* end namespace lama */

} /* end namespace scai */


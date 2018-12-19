/**
 * @file fft.hpp
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
 * @brief Fast fourier transformation for vector and matrices
 * @author Thomas Brandes
 * @date 21.03.2018
 */
#pragma once

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/utilskernel/FFTUtils.hpp>

#include <scai/tracing.hpp>

namespace scai
{

namespace lama
{

#ifdef SCAI_COMPLEX_SUPPORTED

/** 
 *  @brief Apply FFT to a one-dimensional vector
 *
 *  - FFT can only be applied to complex vectors
 *  - The size must be a power of 2
 *  - use resize function to fill up or truncate matrices
 */
template<typename ComplexType>
void fft1D( DenseVector<ComplexType>& denseVector, const int direction )
{
    if (  !denseVector.getDistribution().isReplicated() )
    {
        // FFT not available yet for distributed vectors, so replicate it

        auto repDist = std::make_shared<dmemo::NoDistribution>( denseVector.size() );
        auto saveDist = denseVector.getDistributionPtr();
        denseVector.redistribute( repDist );
        fft1D( denseVector, direction );
        denseVector.redistribute( saveDist );
        return;
    }

    hmemo::ContextPtr ctx = denseVector.getContextPtr();   // preferred location for FFT
    
    hmemo::HArray<ComplexType>& localData = denseVector.getLocalValues();

    const IndexType size = localData.size();

    IndexType m = common::Math::nextpow2( size );
    IndexType n2 = 1 << m;

    SCAI_ASSERT_EQ_ERROR( n2, size, "FFT vector size must be power of 2" )

    const IndexType many = 1;   // one single vector only

    utilskernel::FFTUtils::fftcall<RealType<ComplexType>>( localData, many, size, m, direction, denseVector.getContextPtr() );
}

/** 
 *  @brief Apply forward FFT to a one-dimensional vector
 *
 *  @param[in,out] denseVector vector to which (forward) FFT is applied to, size must be a power of 2
 *
 */
template<typename ComplexType>
void fft( DenseVector<ComplexType>& denseVector )
{
    fft1D( denseVector, 1 );
}

/** 
 *  @brief Apply inverse FFT to a one-dimensional vector
 *
 *  @param[in,out] denseVector vector to which (backward) FFT is applied to, size must be a power of 2
 *
 */
template<typename ComplexType>
void ifft( DenseVector<ComplexType>& denseVector )
{
    fft1D( denseVector, -1 );
}

/** 
 *  @brief Apply FFT to each row for a 2D dense matrix
 *
 *  - FFT can only be applied to complex matrices
 *  - The column size must be a power of 2
 *  - use resize function to fill up or truncate matrices
 */
template<typename ComplexType>
void fftRows( 
    DenseMatrix<ComplexType>& data,
    int direction )
{
    IndexType numColumns = data.getNumColumns();  // size of each vector to which FFT is applied

    IndexType m = common::Math::nextpow2( numColumns );
    IndexType n2 = 1 << m;

    SCAI_ASSERT_EQ_ERROR( n2, numColumns, "not power of 2" )

    if ( !data.getColDistribution().isReplicated() )
    {
        SCAI_REGION( "Matrix.fftRowsRep" )
        dmemo::DistributionPtr repColumns = std::make_shared<dmemo::NoDistribution>( data.getNumColumns() );
        dmemo::DistributionPtr distColumns = data.getColDistributionPtr();
        data.redistribute( data.getRowDistributionPtr(), repColumns );
        fftRows( data, direction );
        data.redistribute( data.getRowDistributionPtr(), distColumns );
    }
    else
    {
        SCAI_REGION( "Matrix.fftRows" )

        // fft  for each row, each processor does it on its local rows

        IndexType nLocalRows = data.getRowDistribution().getLocalSize();

        DenseStorage<ComplexType>& localMatrix = data.getLocalStorage();
        hmemo::HArray<ComplexType>& localData = localMatrix.getData();

        utilskernel::FFTUtils::fftcall<RealType<ComplexType>>( localData, nLocalRows, numColumns, m, direction, data.getContextPtr() );
    }
}

template<typename ComplexType>
void fftCols(
    DenseMatrix<ComplexType>& data,
    int direction )
{
    // transpose the matrix, apply FFT row-wise and transpose back

    if ( !data.getRowDistribution().isReplicated() && data.getColDistribution().isReplicated() )
    {
        // distributed rows, columns must also be distributed for transpose

        dmemo::DistributionPtr saveColDistributionPtr = data.getColDistributionPtr();
        dmemo::DistributionPtr blockDist = std::make_shared<dmemo::BlockDistribution>( data.getNumColumns() );
        data.redistribute( data.getRowDistributionPtr(), blockDist );
        fftCols( data, direction );
        data.redistribute( data.getRowDistributionPtr(), saveColDistributionPtr );
    }
    else if ( data.getRowDistribution().isReplicated() && !data.getColDistribution().isReplicated() )
    {
        // replicated rows, columns must also be replicated

        dmemo::DistributionPtr saveRowDistributionPtr = data.getRowDistributionPtr();
        dmemo::DistributionPtr blockDist = std::make_shared<dmemo::BlockDistribution>( data.getNumRows() );
        data.redistribute( blockDist, data.getColDistributionPtr() );
        fftCols( data, direction );
        data.redistribute( saveRowDistributionPtr, data.getColDistributionPtr() );
    }
    else
    {   // now transpose works fine, we do it in place

        SCAI_REGION( "Matrix.fftCols" )

        data.assignTranspose( data );
        fftRows( data, direction );
        data.assignTranspose( data );
    }
}

/** 
 *  @brief Apply fast fourier transform to a matrix
 *
 *  - For dim = 0 the FFT is applied to each column
 *  - For dim = 1 the FFT is applied to each row 
 *  - The matrix size in the dimension to which the FFT is applied must be a power of 2
 *  - FFT can only be applied to complex matrices
 *  - use resize function to fill up or truncate matrices
 */

/* --------------------------------------------------------------------- */

template<typename ComplexType>
void fft( 
    DenseMatrix<ComplexType>& data,
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
   else
   {
       SCAI_THROWEXCEPTION( common::InvalidArgumentException, "dim = " << dim << ": illegal, must be 0 or 1" )
   }
}

template<typename ComplexType>
void ifft( 
    DenseMatrix<ComplexType>& data,
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
   else
   {
       SCAI_THROWEXCEPTION( common::InvalidArgumentException, "dim = " << dim << ": illegal, must be 0 or 1" )
   }
}

/**
 * Apply Fast Fourier Transform in both dimensions of a matrix.
 *
 * @param[in,out] data is a two-dimensional dense matrix, must be complex and its sizes a power of 2
 */
template<typename ComplexType>
void fft(
    DenseMatrix<ComplexType>& data )
{
    fftRows( data, 1 );   // fast backward FFT alon rows
    fftCols( data, 1 );   // fast backward FFT along columns
}

/**
 * Apply inverse Fast Fourier Transform in both dimensions of a matrix.
 *
 * @param[in,out] data is a two-dimensional dense matrix, must be complex and its sizes a power of 2
 */
template<typename ComplexType>
void ifft(
    DenseMatrix<ComplexType>& data )
{
    // order does not matter, so fftCols before fftRows works also fine

    fftRows( data, -1 );   // fast backward FFT alon rows
    fftCols( data, -1 );   // fast backward FFT along columns
}

#endif

} /* end namespace lama */

} /* end namespace scai */


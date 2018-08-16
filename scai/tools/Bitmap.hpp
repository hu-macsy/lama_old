/**
 * @file Bitmap.hpp
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
 * @brief Class to spy CSR sparse matrix as a bitmap
 * @author Thomas Brandes
 * @date 17.10.2013
 */

#pragma once

#include <scai/lama.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/io/ImageIO.hpp>

using namespace scai;

template<typename ValueType>
class Bitmap
{
public:

    /** Generate an empty bitmap of a certain size.
     *
     *  @param[in] width of the bitmap
     *  @param[in] height of the bitmap
     */

    Bitmap( const IndexType width, const IndexType height ) :

        mPixelData( lama::GridVector<ValueType>( common::Grid2D( width, height ), ValueType( 0 ) ) )
    {
        // set min and max with the zero element

        mMinVal = common::TypeTraits<ValueType>::getMax();   // neutral element for min operation
        mMaxVal = ValueType( 0 );
    }

    ~Bitmap()
    {
    }

    /** this method fills the pixel array with the sparse pattern */

    void drawCSR( const IndexType nRows, const IndexType nCols, const IndexType ia[], const IndexType ja[], const ValueType[] )
    {
        lama::GridWriteAccess<ValueType> wPixel( mPixelData );

        double multRow = double( wPixel.size( 0 ) ) / double( nRows );
        double multCol = double( wPixel.size( 1 ) ) / double( nCols );

        ValueType one = 1;

        mMinVal = 1 / ( one + one ); // must be > 0 and < 1 

        for ( IndexType i = 0; i < nRows; ++i )
        {
            for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
            {
                const IndexType iRow = static_cast<IndexType>( i * multRow );
                const IndexType iCol = static_cast<IndexType>( ja[j] * multCol );
                wPixel( iRow, iCol ) += one;

                mMaxVal = common::Math::max( mMaxVal, wPixel( iRow, iCol ) );
            }
        }
    }

    void drawCSRValues( const IndexType nRows, const IndexType nCols, const IndexType ia[], const IndexType ja[], const ValueType values[] )
    {
        lama::GridWriteAccess<ValueType> wPixel( mPixelData );

        double multRow = double( wPixel.size( 0 ) ) / double( nRows );
        double multCol = double( wPixel.size( 1 ) ) / double( nCols );

        for ( IndexType i = 0; i < nRows; ++i )
        {
            for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
            {
                ValueType matrixValue = values[j];
                matrixValue = common::Math::abs( matrixValue );

                const IndexType iRow = static_cast<IndexType>( i * multRow );
                const IndexType iCol = static_cast<IndexType>( ja[j] * multCol );
                wPixel( iRow, iCol ) = matrixValue;

                if ( matrixValue > common::TypeTraits<ValueType>::eps0() )
                {
                    mMinVal = common::Math::min( mMinVal, matrixValue );
                }
                mMaxVal = common::Math::max( mMaxVal, matrixValue );
            }
        }

        // diagonals drawn at the end, useful if many matrix entries map to one pixel

        for ( IndexType i = 0; i < nRows; ++i )
        {
            for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
            {
                if ( ja[j] == i )
                {
                    ValueType matrixValue = values[j];
                    matrixValue = common::Math::abs( matrixValue );
                    const IndexType iRow = static_cast<IndexType>( i * multRow );
                    const IndexType iCol = static_cast<IndexType>( i * multCol );
                    wPixel( iRow, iCol ) = matrixValue;
                }
            }
        }
    }

    void write( const std::string outFileName )
    {
        SCAI_ASSERT_GT_ERROR( mMinVal, ValueType( 0 ), 
                              "min value must be positive to guarantee background color for zero entries" )

        lama::ImageIO::writeSC( mPixelData, mMinVal, mMaxVal, outFileName );
    }

private:

    lama::GridVector<ValueType> mPixelData;

    ValueType mMinVal;
    ValueType mMaxVal;
};

/**
 * @file Bitmap.hpp
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

class Bitmap
{
public:

    /** Generate an empty bitmap of a certain size.
     *
     *  @param[in] width of the bitmap
     *  @param[in] height of the bitmap
     */

    Bitmap( const IndexType width, const IndexType height ) :

        mPixelData( lama::GridVector<float>( common::Grid2D( width, height ), -10.0 ) )
    {
        // set min and max with the neutral elements

        mMinVal = common::TypeTraits<float>::getMax();
        mMaxVal = common::TypeTraits<float>::getMin();
    }

    ~Bitmap()
    {
    }

    template<typename ValueType>
    void drawCSR( const IndexType nRows, const IndexType nCols, const IndexType ia[], const IndexType ja[], const ValueType values[] )
    {
        lama::GridWriteAccess<float> wPixel( mPixelData );

        double multRow = double( wPixel.size( 0 ) ) / double( nRows );
        double multCol = double( wPixel.size( 1 ) ) / double( nCols );

        for ( IndexType i = 0; i < nRows; ++i )
        {
            for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
            {
                float matrixValue = static_cast<float>( values[j] );

                const IndexType iRow = static_cast<IndexType>( i * multRow );
                const IndexType iCol = static_cast<IndexType>( ja[j] * multCol );
                wPixel( iRow, iCol ) = matrixValue;

                mMinVal = common::Math::min( mMinVal, matrixValue );
                mMaxVal = common::Math::max( mMaxVal, matrixValue );
            }
        }

        // diagonals drawn at the end

        for ( IndexType i = 0; i < nRows; ++i )
        {
            for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
            {
                if ( ja[j] == i )
                {
                   const IndexType iRow = static_cast<IndexType>( i * multRow );
                   const IndexType iCol = static_cast<IndexType>( i * multCol );
                   wPixel( iRow, iCol ) = static_cast<float>( values[j] );
                }
            }
        }
    }

    void write( const std::string outFileName )
    {
        lama::ImageIO::writeSC( mPixelData, mMinVal, mMaxVal, outFileName );
    }

private:

    lama::GridVector<float> mPixelData;

    float mMinVal;
    float mMaxVal;
};

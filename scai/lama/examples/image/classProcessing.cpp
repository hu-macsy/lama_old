/**
 * @file classProcessing.cpp
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
 * @brief Example program to classify image data
 * @author Thomas Brandes, Dustin Feld
 * @date 05.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

// _Matrix & vector related includes

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc != 3 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <inputFileName> <outputFileName>" << std::endl;
        return -1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    // read in the image file, must be a png/bmp file

    GridVector<ValueType> image;   // size will be ( width , height, ncolors )

    std::cout << "read image as grid vector : " << image << std::endl;

    const common::Grid& grid = image.globalGrid();

    SCAI_ASSERT_EQ_ERROR( 3, grid.nDims(), "not an image" );

    const IndexType height = grid.size( 0 );
    const IndexType width  = grid.size( 1 );

    SCAI_ASSERT_EQ_ERROR( 3, grid.size( 2 ), "not RGB pixels" );

    GridVector<ValueType> dist ( common::Grid2D( height, width ), 0.0f );

    // define the color to be used for classification 

    const ValueType classColor[][3] = { { 0, 128, 128 }, { 128, 0, 128 }, { 128, 128, 0 },
                                    { 0, 128, 0 }, { 128, 0, 0}, { 0, 0, 128 } };

    IndexType nClassColors = sizeof( classColor ) / ( 3 * sizeof( ValueType ) );

    GridVector<ValueType> classImage( grid, 0 );

    // Idea: set for each pixel the classification color that is the closest

    std::cout << "Class image = " << classImage << std::endl;

    {
        GridWriteAccess<ValueType> wClass( classImage );
        GridReadAccess<ValueType> rImage( image );
        GridWriteAccess<ValueType> wDist( dist );

        for ( IndexType i = 0; i < height; i++ )
        {
            for ( IndexType j = 0; j < width; j++ )
            {
               ValueType dist = 0.0;  // distance of image(i,j) to classColor[0] 

               for ( IndexType c = 0; c < 3; ++c )
               {
                   ValueType diff = classColor[0][c] - rImage( i, j, c );
                   dist += diff * diff;
               }

               dist = common::Math::sqrt( dist );

               wDist( i, j ) = dist;

               wClass( i, j, 0 ) = classColor[0][0];
               wClass( i, j, 1 ) = classColor[0][1];
               wClass( i, j, 2 ) = classColor[0][2];
            }
        }
    }

    /* Pseudo code for next block:

       colDist( i, j ) = l2norm( image(i,j,:) - col(k) 
       where ( colDist(:,:) < dist(:,:)
           classImage = col[k]
           dist  = colDist
        endwhere
    */

    for ( IndexType k = 1; k < nClassColors; ++k )
    {
        GridWriteAccess<ValueType> wClass( classImage );
        GridReadAccess<ValueType> rImage( image );
        GridWriteAccess<ValueType> wDist( dist );

        for ( IndexType i = 0; i < height; i++ )
        {
            for ( IndexType j = 0; j < width; j++ )
            {
               ValueType dist = 0.0;

               for ( IndexType c = 0; c < 3; ++c )
               {
                   ValueType diff = classColor[k][c] - rImage( i, j, c );
                   dist += diff * diff;
               }

               dist = common::Math::sqrt( dist );

               if ( dist < wDist( i, j ) )
               {
                   wDist( i, j ) = dist;

                   wClass( i, j, 0 ) = classColor[k][0];
                   wClass( i, j, 1 ) = classColor[k][1];
                   wClass( i, j, 2 ) = classColor[k][2];
               }
            }
        }
    }

    if ( false )
    {
        // apply edge detection

        ValueType edge_pattern[9] = { 0, -1, 0, -1, 4, -1, 0, -1, 0 };

        common::Stencil3D<ValueType> stencil( 3, 3, 1, edge_pattern );

        StencilMatrix<ValueType> edgeIt( image.getDistributionPtr(), stencil );

        GridVector<ValueType> imageNew;

        imageNew = edgeIt * classImage;

        classImage.swap( imageNew );
    }

    std::cout << "new image = " << classImage << std::endl;

    classImage.writeToFile( outputFileName );
}

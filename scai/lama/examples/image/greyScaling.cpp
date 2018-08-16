/**
 * @file greyScaling.cpp
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
 * @brief Example program to greyscale an image file
 * @author Thomas Brandes
 * @date 10.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridSection.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

// _Matrix & vector related includes

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

/** These are the scale factor to combine the RGB values */

const ValueType greyWeights[] = {  0.2989, 0.5870, 0.1140 };

/** Computing a greyscaled image via grid read/write accesses.
 *
 *  This solution uses the classes GridReadAccess and GridWriteAccess
 *  that allow a very convenient multi-dimensional addressing of the local part
 *  of the grid vector.
 */

template<typename ValueType> 
void greyScale1( GridVector<ValueType>& imageGrey, const GridVector<ValueType>& image )
{
    const common::Grid& grid = image.globalGrid();

    const IndexType height = grid.size( 0 );
    const IndexType width  = grid.size( 1 );

    // output image gets same size, shape, distribution

    imageGrey.allocate( image.getDistributionPtr() );

    {
        GridWriteAccess<ValueType> wImageGrey( imageGrey );
        GridReadAccess<ValueType> rImage( image );

        for ( IndexType i = 0; i < height; ++i )
        {
            for ( IndexType j = 0; j < width; ++j )
            {
                ValueType g = rImage( i, j, 0 ) * greyWeights[0] +
                              rImage( i, j, 1 ) * greyWeights[1] +
                              rImage( i, j, 2 ) * greyWeights[2];

                wImageGrey( i, j, 0 ) = g;
                wImageGrey( i, j, 1 ) = g;
                wImageGrey( i, j, 2 ) = g;
            }
        }
    }
}

template<typename ValueType> 
void greyScale2( GridVector<ValueType>& imageGrey, const GridVector<ValueType>& image )
{
    // grey scaling can be formulated as a stencil operation in the 3-rd dimension

    common::Stencil3D<ValueType> greyStencil( 1, 1, 3, greyWeights );

    StencilMatrix<ValueType> greyMatrix( image.getDistributionPtr(), greyStencil );

    imageGrey = greyMatrix * image;

    // define sections for the different colors R, G, and B

    GridSection<ValueType>sectionR ( imageGrey, Range(), Range(), 0 );
    GridSection<ValueType>sectionG ( imageGrey, Range(), Range(), 1 );
    GridSection<ValueType>sectionB ( imageGrey, Range(), Range(), 2 );
    
    // due to the stencil only the G of RGB values contains the right value

    sectionR = sectionG;
    sectionB = sectionG;

}

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

    // read in the image file, must be a png file

    GridVector<ValueType> image( inputFileName );   // size will be ( width , height, ncolors )

    std::cout << "read image as grid vector : " << image << std::endl;

    const common::Grid& grid = image.globalGrid();

    SCAI_ASSERT_EQ_ERROR( 3, grid.nDims(), "not an image" );
    SCAI_ASSERT_EQ_ERROR( 3, grid.size( 2 ), "not RGB pixels" );

    GridVector<ValueType> imageGrey;

    // greyScale1( imageGrey, image );
    greyScale2( imageGrey, image );
 
    imageGrey.writeToFile( outputFileName );
}

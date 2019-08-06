/**
 * @file imageProcessing.cpp
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
 * @brief Example program to work on image data
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/common/Stencil.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

#include "ConvolutionMatrices.hpp" 

// _Matrix & vector related includes

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

void gaussianBlur(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage, const IndexType radius, const ValueType sigma){

    // low pass filter
    // radius = amount of blur, sigma = standard deviation

    common::Stencil3D<ValueType> stencil;
    
    ValueType num = 0;
    for(IndexType i = -radius; i<radius+1; i++)
    {
        for(IndexType j = -radius; j<radius+1; j++)
        {
            num = num +(1/(2*M_PI*(sigma*sigma))*pow(M_E, -((i*i+j*j)/(2*sigma*sigma))));
        }

    }

    // use num so that sum(all points of stencil)=1

    for(IndexType i = -radius; i<radius+1; i++)
    {
        for(IndexType j = -radius; j<radius+1; j++)
        {
            stencil.addPoint(i,j,0,((1/(2*M_PI*(sigma*sigma))*pow(M_E, -((i*i+j*j)/(2*sigma*sigma)))))/num);
        }

    }

    StencilMatrix<ValueType> blur( inputImage.getDistributionPtr(), stencil);

    outputImage = blur*inputImage;
}

void meanFilter(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage, IndexType radius){

    // low pass filter
    // radius = amount of blur

    ValueType mean =(radius+radius+1)*(radius+radius+1);
    mean = 1/mean; 

    common::Stencil3D<ValueType> stencil;

    for(IndexType i = -radius; i<radius+1; i++)
    {
        for(IndexType j = -radius; j<radius+1; j++)
        {
            stencil.addPoint(i,j,0,mean);
        }

    }

    StencilMatrix<ValueType> matrix( inputImage.getDistributionPtr(), stencil );

    outputImage = matrix* inputImage;

}


void grayScale(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage){

    const ValueType red   = 0.3;
    const ValueType green = 0.59;
    const ValueType blue  = 0.11;
    const ValueType sumOfWeights = red + green + blue;

    SCAI_ASSERT_EQ_ERROR( sumOfWeights, 1, "redWeighting + greenWeighting + blueWeighting = " << sumOfWeights << " must be 1" )

    outputImage.allocate( inputImage.getDistributionPtr() );
    const common::Grid& grid = inputImage.globalGrid();

    GridWriteAccess<ValueType> wImage( outputImage );
    GridReadAccess<ValueType>  rImage( inputImage );

    for(IndexType i = 0 ; i<grid.size(0); i++)
    {
        for(IndexType j = 0 ; j<grid.size(1); j++)
        {
            ValueType mean =( red * rImage(i,j,0) + green*rImage(i,j,1) + blue*rImage(i,j,3) );

            wImage( i,  j,  0) = mean;
            wImage( i,  j,  1) = mean;
            wImage( i,  j,  2) = mean;
        }
    } 
}

void cutOffColour(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage){

    // the range of values of RGB is [0,255], cut off the higher/lower Values

    const common::Grid& grid = inputImage.globalGrid();
    outputImage=inputImage;

    GridWriteAccess<ValueType> wImage( outputImage );
    GridReadAccess<ValueType>  rImage( inputImage );

    for(IndexType i = 0 ; i<grid.size(0); i++) 
    {
        for(IndexType j = 0 ; j<grid.size(1); j++)
        {
            for(IndexType x = 0 ; x<3 ; x++)
            {
                if (rImage(i,j,x)<0)
                {
                    wImage(i,j,x)=0;
                }
                else if (rImage(i,j,x)>255)
                {
                    wImage(i,j,x)=255;
                }
            }
        }
    }

}

void sobelFilter(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage){

    // sobel filter can be used for Edge detection

    common::Stencil3D<ValueType> stencilSobelX ( 3, 3, 1, imageprocessing::sobelX );
    common::Stencil3D<ValueType> stencilSobelY ( 3, 3, 1, imageprocessing::sobelY );

    StencilMatrix<ValueType> matrixSobelX( inputImage.getDistributionPtr(), stencilSobelX );
    StencilMatrix<ValueType> matrixSobelY( inputImage.getDistributionPtr(), stencilSobelY );

    GridVector<ValueType> SobelX;
    GridVector<ValueType> SobelY;
    GridVector<ValueType> Gray;

    SobelX  =  matrixSobelX * inputImage;
    SobelY  =  matrixSobelY * inputImage;
    SobelX.dotProduct(SobelY); // both x and y filters are applied

    grayScale(Gray, SobelX);

    cutOffColour(outputImage, Gray);

}

void unsharpMask(GridVector<ValueType>& outputImage, const GridVector<ValueType>& inputImage, const ValueType amount, const IndexType radius){

    // the unsharp mask is used to sharpen an image by using a blurred copy
    // amount = strength of sharpening
    // radius = amount of blur

    const ValueType sigma = 1.5;

    GridVector<ValueType> blur;
    GridVector<ValueType> mask;
    GridVector<ValueType> shapen;


    mask.allocate( inputImage.getDistributionPtr() );

    gaussianBlur(blur, inputImage, radius, sigma);
    mask= inputImage-blur;
    cutOffColour(shapen, mask);

    shapen = shapen * amount + inputImage;
    cutOffColour(outputImage, shapen);
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

    auto image = read<GridVector<ValueType>>( inputFileName ); 

    std::cout << "read image as grid vector : " << image << std::endl;

    SCAI_ASSERT_EQ_ERROR( image.nDims(), 3, "no color image data" );

    const common::Grid& grid = image.globalGrid(); 
    SCAI_ASSERT_EQ_ERROR( 3, grid.size( 2 ), "not RGB pixels" );

    // apply stencil on the pixels, do not apply on the colors in 3-rd dimension 

    common::Stencil3D<ValueType> stencil( 5, 5, 1, imageprocessing::findEdges );

    // common::Stencil3D<ValueType> stencil( 3, 3, 1, imageprocessing::blur );


    // common::Stencil3D<ValueType> stencil( 3, 3, 1, imageprocessing::sharpen );


    // common::Stencil3D<ValueType> stencil( 1 );   // one-point stencil just as dummy


    StencilMatrix<ValueType> m( image.getDistributionPtr(), stencil );

    std::cout << "stencil matrix = " << m << std::endl;

    GridVector<ValueType> imageNew; 

    imageNew  =  m * image;

    //gaussianBlur(imageNew, image, 6, 1.5); 
    //grayScale(imageNew,image);
    //cutOffColour(imageNew,image);
    //sobelFilter(imageNew, image);
    //meanFilter(imageNew,image,6);
    //unsharpMask(imageNew,image,0.5,4);

    std::cout << "new image = " << imageNew << std::endl;


    imageNew.writeToFile( outputFileName );
}

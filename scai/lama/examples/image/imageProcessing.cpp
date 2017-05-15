/**
 * @file imageProcessing.cpp
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

// Matrix & vector related includes

using namespace scai;
using namespace lama;


GridVector<float>  gaussianBlur(GridVector<float> inputImage, IndexType radius, float sigma){

    // radius = amount of blur, sigma = standard deviation

    common::Stencil3D<float> stencil;
    
    float num = 0;
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

    StencilMatrix<float> blur( inputImage.getDistributionPtr(), stencil);

    GridVector<float> output;
    output = blur*inputImage;

    return output;

}

GridVector<float> grayScale(const GridVector<float>& inputImage){

    const float red   = 0.3;
    const float green = 0.59;
    const float blue  = 0.11;
    const float sumOfWeights = red + green + blue;

    if(sumOfWeights != 1)
    {
        SCAI_ASSERT_EQ_ERROR( sumOfWeights, 1, "redWeighting + greenWeighting + blueWeighting = " << sumOfWeights << " must be 1" )
    }

    const common::Grid& grid = inputImage.globalGrid();

    GridVector<float> outputImage;
    outputImage.allocate( inputImage.getDistributionPtr() );

    GridWriteAccess<float> wImage( outputImage );
    GridReadAccess<float>  rImage( inputImage );

    for(IndexType i = 0 ; i<grid.size(0); i++)
    {
        for(IndexType j = 0 ; j<grid.size(1); j++)
        {
            Scalar mean =( red*rImage(i,j,0) + green*rImage(i,j,1) + blue*rImage(i,j,3) );

            wImage( i,  j,  0) = mean.getValue<float>();
            wImage( i,  j,  1) = mean.getValue<float>();
            wImage( i,  j,  2) = mean.getValue<float>();
        }
    } 
    return outputImage;
}

GridVector<float> cutOffColour(const GridVector<float>& inputImage){

    // the range of values of RGB is [0,255], cut off the higher/lower Values

    const common::Grid& grid = inputImage.globalGrid();

    GridVector<float> outputImage;
    outputImage.allocate( inputImage.getDistributionPtr() );

    GridWriteAccess<float> wImage( outputImage );
    GridReadAccess<float>  rImage( inputImage );

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

    return outputImage;
}

GridVector<float>  sobelFilter(GridVector<float> inputImage){

    // sobel filter can be used for Edge detection

    common::Stencil3D<float> stencilSobelX ( 3, 3, 1, imageprocessing::sobelX );
    common::Stencil3D<float> stencilSobelY ( 3, 3, 1, imageprocessing::sobelY );

    StencilMatrix<float> matrixSobelX( inputImage.getDistributionPtr(), stencilSobelX );
    StencilMatrix<float> matrixSobelY( inputImage.getDistributionPtr(), stencilSobelY );

    GridVector<float> SobelX;
    GridVector<float> SobelY;
    GridVector<float> Gray;

    SobelX  =  matrixSobelX * inputImage;
    SobelY  =  matrixSobelY * inputImage;
    SobelX.dotProduct(SobelY); // both x and y filters are applied

    Gray = grayScale(SobelX);

    GridVector<float> outputImage= cutOffColour(Gray);

    return outputImage;
}

GridVector<float>  meanFilter(GridVector<float> inputImage, IndexType radius){
    // radius = = amount of blur

    float mean =(radius+radius+1)*(radius+radius+1);
    mean = 1/mean; 

    common::Stencil3D<float> stencil;

    for(IndexType i = -radius; i<radius+1; i++)
    {
        for(IndexType j = -radius; j<radius+1; j++)
        {
            stencil.addPoint(i,j,0,mean);
        }

    }

    StencilMatrix<float> matrix( inputImage.getDistributionPtr(), stencil );

    GridVector<float> outputImage;
    outputImage.allocate( inputImage.getDistributionPtr() );
    outputImage = inputImage * matrix;

    return outputImage;
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

    GridVector<float> image;   // size will be ( width , height, ncolors )

    ImageIO::read( image, inputFileName );

    std::cout << "read image as grid vector : " << image << std::endl;

    SCAI_ASSERT_EQ_ERROR( image.nDims(), 3, "no color image data" )

    // apply stencil on the pixels, do not apply on the colors in 3-rd dimension 

     common::Stencil3D<float> stencil( 5, 5, 1, imageprocessing::findEdges );

    // common::Stencil3D<float> stencil( 3, 3, 1, imageprocessing::blur );


    // common::Stencil3D<float> stencil( 3, 3, 1, imageprocessing::sharpen );


    // common::Stencil3D<float> stencil( 1 );   // one-point stencil just as dummy


    StencilMatrix<float> m( image.getDistributionPtr(), stencil );

    std::cout << "stencil matrix = " << m << std::endl;

    GridVector<float> imageNew;

    imageNew  =  m * image;

    //GridVector<float> imageGausBlur= gaussianBlur(image, 6, 1.5);
    //GridVector<float> imageGrey= grayScale(image);
    //GridVector<float> imagesobelFilter= sobelFilter(image);
    //GridVector<float> imageblur = meanFilter(image,1);
    // std::cout << "new image = " << imageNew << std::endl;


    ImageIO::write( imageNew, outputFileName );
}

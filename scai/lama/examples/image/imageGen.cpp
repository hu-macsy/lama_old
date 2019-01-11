/**
 * @file imageGen.cpp
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
 * @brief Example program to demonstrate image generation and indexing via GridVector and GridSection.
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/common/Stencil.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/lama/GridSection.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>


using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc < 2 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <outputFileName> [ width [ height ] ]" << std::endl;
        return -1;
    }

    std::string outputFileName = argv[1];

    IndexType width =  150;
    IndexType height = 100;

    if ( argc > 2 )
    {
        std::istringstream input( argv[2] );
        input >> width;
        height = width;

        if ( argc > 3 )
        {
            std::istringstream input( argv[3] );
            input >> height;
        }
    }

    IndexType border = common::Math::min( width, height ) / 2 - 5;

    // A three-dimensional grid with 3 entries in the last dimension can be
    // used as an image.
    //
    //  ATTENTION:   specify height x width  ( conform with nrows x ncols )
    //

    GridVector<ValueType> image( common::Grid3D( height, width, 3 ) );

    // Indexing of image is like indexing a matrix + entry for color
    //
    //   ( 0, 0 )  ( 0, 1 )  ( 0, 2 ) ---->
    //   ( 1, 0 )  ( 1, 1 )  ( 1, 2 ) ---->
    //   ( 2, 0 )  ( 2, 1 )  ( 2, 2 ) ---->
    //       |         |         |

    // top left : white
    image( Range( 0, border ), Range( 0, border ), 0 ) = 255;
    image( Range( 0, border ), Range( 0, border ), 1 ) = 255;
    image( Range( 0, border ), Range( 0, border ), 2 ) = 255;

    // top right : blue
    image( Range( 0, border ), Range( width - border, width ), 2 ) = 255;

    // bottom left: red
    image( Range( height - border, height), Range( 0, border ), 0 ) = 255;

    // bottom left: green
    image( Range( height - border, height), Range( width - border, width ), 1 ) = 255;

    // ImageIO::write( image, outputFileName );

    try 
    {
        image.writeToFile( outputFileName );
        std::cout << "Image successfully written to file " << outputFileName << std::endl;
    }
    catch ( common::Exception& ex )
    {
        std::cout << "Failed to write image: " << ex.what() << std::endl;
    }
}

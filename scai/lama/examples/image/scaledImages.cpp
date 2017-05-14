/**
 * @file scaledImages.cpp
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
 * @brief Example program to generate scaled images
 * @author Thomas Brandes
 * @date 14.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/common/Settings.hpp>

// Matrix & vector related includes

using namespace scai;
using namespace lama;


int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc != 1 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << std::endl;
        return -1;
    }

    const IndexType M = 1920;
    const IndexType N = 1024;

    GridVector<float> arrayData( common::Grid2D( M, N ), 0 );

    {
        GridWriteAccess<float> wArray( arrayData );

        for ( IndexType i = 0; i < M; ++i ) 
        {
            for ( IndexType j = 0; j < N; ++j )
            {
                wArray( i, j ) = static_cast<float>( M + N );
            }
        }
    }

    ImageIO::writeSC( arrayData, "array.png" );
}

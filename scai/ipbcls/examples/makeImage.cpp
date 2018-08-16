/**
 * @file makeImage.cpp
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
 * @brief Convert vector data to a scaled image
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#include <scai/common/Settings.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/io/ImageIO.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc < 5 )
    {
        std::cout << "Please call " << argv[0] 
                  << " <vector_filename> nx ny <image_filename>" << std::endl;
        return -1;
    }

    IndexType nx = atoi( argv[2] );
    IndexType ny = atoi( argv[3] );

    auto x = read<DenseVector<ValueType>>( argv[1] );

    SCAI_ASSERT_EQ_ERROR( nx * ny, x.size(), "vector file does not match size " << nx << " x " << ny )

    IndexType scale = 4;

    GridVector<ValueType> xData( common::Grid2D( scale * nx, scale * ny ), 0 );

    {
        hmemo::ReadAccess<ValueType> rX( x.getLocalValues() );
        GridWriteAccess<ValueType> wXData( xData );

        for ( IndexType i = 0; i < nx; ++i )
        {
            for ( IndexType j = 0; j < ny; ++j )
            {
                for ( IndexType i1 = 0; i1 < scale; ++i1 )
                {
                    for ( IndexType j1 = 0; j1 < scale; ++j1 )
                    {
                        wXData( i * scale + i1 , j * scale + j1 ) = rX[ i * ny + j ];
                    }
                }
            }
        }
    }

    ImageIO::writeSC( xData, argv[4] );

    std::cout << "Written autoscale image " << nx << " x " << ny << " to file " << argv[4] << std::endl;
}

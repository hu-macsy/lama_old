/**
 * @file makeImage.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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

/**
 * @file gridIO.cpp
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
 * @date 19.05.2017
 */

#include <scai/lama/io/ImageIO.hpp>

#include <scai/common/Stencil.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;

/**
 *  @brief Example program that writes a grid vector to a file and reads it back from the file
 */
int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments:
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 2 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <fileName>" << std::endl;
        return -1;
    }

    std::string fileName = argv[1];

    const IndexType n1 = 5;
    const IndexType n2 = 3;
    const IndexType n3 = 2;
    const IndexType n4 = 4;

    auto gv1 = GridVector<double>( common::Grid4D( n1, n2, n3, n4 ), 0 );

    {
        GridWriteAccess<double> wGV1( gv1 );

        for ( IndexType i1 = 0; i1 < n1; ++i1 )
        {
            for ( IndexType i2 = 0; i2 < n2; ++i2 )
            {
                for ( IndexType i3 = 0; i3 < n3; ++i3 )
                {
                    for ( IndexType i4 = 0; i4 < n4; ++i4 )
                    {
                        wGV1( i1, i2, i3, i4 ) = 1000 * ( i1 + 1 ) + 100 * ( i2 + 1 ) + 10 * ( i3 + 1 ) + i4 + 1;
                    }
                }
            }
        }
    }

    std::cout << "write to file " << fileName << " this grid vector: " << gv1 << std::endl;

    gv1.writeToFile( fileName );

    GridVector<double> gv2 = read<GridVector<double>>( fileName );

    std::cout << "read from file " << fileName << " this grid vector: " << gv2 << std::endl;

    gv2.replicate();

    SCAI_ASSERT_EQ_ERROR( gv1.globalGrid(), gv2.globalGrid(), "mismatch grid size" );

    {
        GridReadAccess<double> rGV1( gv1 );
        GridReadAccess<double> rGV2( gv2 );

        for ( IndexType i = 0; i < gv1.size(); ++i )
        {
            SCAI_ASSERT_EQ_ERROR( rGV1[i], rGV2[i], "different val at i = " << i )
        }
    }
}

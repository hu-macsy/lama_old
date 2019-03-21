/**
 * @file denseIO.cpp
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
 * @brief Example program to read/write dense matrix as 2D grid vector
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

    IndexType m = 7;
    IndexType n = 5;

    DenseMatrix<double> dMatrix( m, n );

    {
        hmemo::WriteAccess<double> wDense( dMatrix.getLocalStorage().getData() );

        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = 0; j < n; ++j )
            {
                wDense[ i * n + j ] = static_cast<double>( 10 * i + j ) ;
            }
        }
    }

    std::cout << dMatrix << ": write it to file " << fileName << std::endl;

    dMatrix.writeToFile( fileName );

    std::cout << dMatrix << " written to file " << fileName << std::endl;

    auto dg = read<GridVector<double>>( fileName );

    std::cout << "read this grid vector from file " << fileName << ": " << dg << std::endl;

    SCAI_ASSERT_EQ_ERROR( dg.globalGrid(), common::Grid2D( m, n ), "mismatch" );

    // dg might be distributed, for check we replicate it

    dg.replicate();

    {
        hmemo::ReadAccess<double> rGrid( dg.getLocalValues() );
        hmemo::ReadAccess<double> rDense( dMatrix.getLocalStorage().getData() );

        for ( IndexType i = 0; i < dg.size(); ++i )
        {
            SCAI_ASSERT_EQ_ERROR( rGrid[i], rDense[i], "different val at i = " << i )
        }
    }
}

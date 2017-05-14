/**
 * @file gridExample.cpp
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
 * @brief Example program to work on grid data using reductions
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;

using namespace common;

typedef double ValueType;
typedef ComplexDouble ComplexType;

SCAI_LOG_DEF_LOGGER( logger, "main" )

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments:
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 1 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <inputFileName> <outputFileName>" << std::endl;
        return -1;
    }

    const IndexType n = 4;
    const IndexType m = 5;
    const IndexType p = 6;

    GridVector<double> A( Grid2D( m, p ), 1.0 );
    GridVector<double> B( Grid2D( p, n ), 1.0 );

    {
        GridWriteAccess<double> wA( A );

        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType k = 0; k < p; ++k )
            {
                wA( i, k ) = 10 * i - k;
            }
        }
    }

    {
        GridWriteAccess<double> wB( B );

        for ( IndexType k = 0; k < p; ++k )
        {
            for ( IndexType j = 0; j < n; ++j )
            {
                wB( k, j ) = 2 * k + 3 * j;
            }
        }
    }

    GridVector<double> C( Grid2D( m, n ), 0.0 );
    GridVector<double> C1( Grid2D( m, n ), 0.0 );

    C.gemm( 1, A, B );

    {
        GridReadAccess<double> rA( A );
        GridReadAccess<double> rB( B );
        GridWriteAccess<double> wC( C1 );

        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = 0; j < n; ++j )
            {
                for ( IndexType k = 0; k < p; ++k )
                {
                    wC( i, j ) += rA( i, k ) * rB( k, j );
                }
            }
        }
    }

    C.writeToFile( "C.mtx" );
    C1.writeToFile( "C1.mtx" );

    std::cout << "diff = " << C.maxDiffNorm( C1 ) << std::endl;
}

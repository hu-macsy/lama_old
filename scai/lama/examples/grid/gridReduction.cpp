/**
 * @file gridReduction.cpp
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

typedef double ValueType;

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

    const IndexType nx   = 5;
    const IndexType nz   = 5;
    const IndexType nsrc = 10;
    const IndexType nf   = 5;

    GridVector<ValueType> Pplus( common::Grid4D( nx, nz, nsrc, nf ), 1 );

    // %% SUMMATION recalc r
    // R=sum(sum(abs(Pplus).^2,4),3);

    GridVector<ValueType> R1( common::Grid3D( nx, nz, nsrc ), 0 );
    GridVector<ValueType> R( common::Grid2D( nx, nz ), 0 );
    GridVector<ValueType> tmp;

    tmp = R * R;

    // ToDo:

    // R1.reduce( tmp, 3, common::binary::ADD );
    // R.reduce( R1, 2, common::binary::ADD );

    {
        GridWriteAccess<ValueType> wR( R );
        GridReadAccess<ValueType> rPplus( Pplus );

        for ( IndexType iz = 0; iz < nz; ++iz )
        {
            for ( IndexType ix = 0; ix < nx; ++ix )
            {
                ValueType sum = 0;

                for ( IndexType isrc = 0; isrc < nsrc; ++isrc )
                {
                    for ( IndexType i_f = 0; i_f < nf; ++i_f )
                    {
                        ValueType s = common::Math::abs( rPplus( iz, ix, isrc, i_f ) );
                        sum = sum + s * s;
                    }
                }
    
                wR( ix, iz ) = sum;
            }
        }
    }
}

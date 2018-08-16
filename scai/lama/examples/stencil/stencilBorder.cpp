/**
 * @file lama/examples/stencil/stencilBorder.cpp
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
 * @brief Demo of stencil class and generation of Stencil matrix
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/NoCommunicator.hpp>

// import common 
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace dmemo;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    common::Stencil1D<ValueType> stencilFD8;

    stencilFD8.reserve( 8 );   // just for convenience, not mandatory

    stencilFD8.addPoint( -3, -11 );
    stencilFD8.addPoint( -2, +7 );
    stencilFD8.addPoint( -1, -3 );
    stencilFD8.addPoint( 0, 1 );
    stencilFD8.addPoint( 1, -1 );
    stencilFD8.addPoint( 2, 3 ) ;
    stencilFD8.addPoint( 3, -7 );
    stencilFD8.addPoint( 4, 11 );

    // 1-dimensional stencils can be combined

    const IndexType N1 = 10;

    common::Grid1D grid( N1 );

    grid.setBorderType( 0, common::Grid::BORDER_ABSORBING, common::Grid::BORDER_ABSORBING );

    StencilStorage<ValueType> st( grid, stencilFD8 );

    for ( IndexType i = 0; i < N1; ++i )
    {
        HArray<ValueType> vals;
        HArray<IndexType> ja;

        st.getSparseRow( ja, vals, i );

        SCAI_ASSERT_EQ_ERROR( ja.size(), vals.size(), "serious mismatch"  );

        hmemo::ReadAccess<ValueType> rValues( vals );
        hmemo::ReadAccess<IndexType> rJA( ja );

        std::cout << "Row " << i << " : ";

        for ( IndexType jj = 0; jj < ja.size(); ++jj )
        {
            std::cout << " " << rJA[jj] << ":" << rValues[jj];
        }

        std::cout << std::endl;
    }

    auto csr = convert<CSRStorage<ValueType>>( st );

    for ( IndexType i = 0; i < N1; ++i )
    {
        HArray<ValueType> vals;
        HArray<IndexType> ja;

        csr.getSparseRow( ja, vals, i );

        SCAI_ASSERT_EQ_ERROR( ja.size(), vals.size(), "serious mismatch"  );

        hmemo::ReadAccess<ValueType> rValues( vals );
        hmemo::ReadAccess<IndexType> rJA( ja );

        std::cout << "Row " << i << " : ";

        for ( IndexType jj = 0; jj < ja.size(); ++jj )
        {
            std::cout << " " << rJA[jj] << ":" << rValues[jj];
        }

        std::cout << std::endl;
    }
}

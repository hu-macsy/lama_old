/**
 * @file kaczmarz.cpp
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
 * @brief Kaczmarz method with sparse vector
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::hmemo;
using namespace scai::lama;
using namespace scai::dmemo;

int main()
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    SCAI_REGION( "main.Kacmarz" )

    CSRSparseMatrix<double> sMatrix( "matrix.frm" );
    DenseMatrix<double> dMatrix( "matrix.frm" );

    IndexType size = sMatrix.getNumRows();

    DistributionPtr blockDist( new BlockDistribution( size, comm ) );

    sMatrix.redistribute( blockDist, blockDist );
    dMatrix.redistribute( blockDist, blockDist );

    std::cout << *comm << ": Matrix: " << sMatrix << ", " << dMatrix << std::endl;

    SparseVector<double> sRow;
    DenseVector<double> dRow;

    for ( IndexType i = 0; i < size; ++i )
    {
        sMatrix.getRow( sRow, i );
        dMatrix.getRow( dRow, i );

        Scalar s = sRow.dotProduct( sRow );
        Scalar d = dRow.dotProduct( dRow );

        SCAI_ASSERT_EQ_ERROR( s, d, "Error iter i = " << i ) 
    }
}

/**
 * @file example.cpp
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
 * @brief Example program used for LAMA module in Users Guide.
 * @author The LAMA development team
 * @date 17.05.2014
 */

#include <scai/lama.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

int main()
{
   auto csrMatrix = read<CSRSparseMatrix<ValueType>>( "gr_30_30.mtx" );

   SCAI_ASSERT_EQ_ERROR( csrMatrix.getNumRows(), csrMatrix.getNumColumns(), "input matrix not square" )

   const IndexType size = csrMatrix.getNumRows();

   auto dist = dmemo::blockDistribution( size );

   csrMatrix.redistribute( dist, dist );

   auto vector = denseVector<ValueType>( dist, 1 );
   auto result = denseVectorEval( csrMatrix * vector );

   result.writeToFile( "result.mtx" );
}

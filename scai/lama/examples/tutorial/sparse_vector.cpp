/**
 * @file sparse_vector.cpp
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
 * @brief Demo program of how to use SparseVector efficiently
 * @author 
 * @date 17.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>

#include <scai/dmemo/NoDistribution.hpp>

#include <iostream>

using namespace scai;
using namespace lama;
using hmemo::HArray;

int main()
{
    auto ctx = hmemo::Context::getContextPtr();

    /** Take default real type for this example. */

    typedef DefaultReal ValueType;

    const IndexType N = 1000;
    const IndexType NIter = 50;

    auto dist = dmemo::noDistribution( N );

    auto vec = denseVectorLinear<ValueType>( dist, 1, 1 );

    hmemo::HArray<IndexType> sparseIndexes( { 3, 8, 11, 19, 20, N - 20, N - 10, N -1 }, ctx );
    hmemo::HArray<ValueType> sparseValues1( { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 }, ctx );
    hmemo::HArray<ValueType> sparseValues2( { 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08 }, ctx );
    hmemo::HArray<ValueType> sparseInv( { 0.8, 0.9, 0.95, 0.99, 1.01, 1.05, 1.1, 1.2 }, ctx );

    const IndexType nnz = sparseIndexes.size();

    SparseVector<ValueType> a( dist, sparseIndexes, sparseValues1, 0, ctx );
    SparseVector<ValueType> b( dist, sparseIndexes, sparseValues2, 0, ctx );
    SparseVector<ValueType> psi( dist, sparseIndexes, HArray<ValueType>( nnz, ValueType( 0 ), ctx ), 0, ctx );
    SparseVector<ValueType> inv( dist, sparseIndexes, HArray<ValueType>( nnz, ValueType( 1 ), ctx ), 1, ctx );
   
    SparseVector<ValueType> tmp( ctx );

    for ( IndexType i = 0; i < NIter; ++i )
    {
        tmp = a;
        psi *= b; 
        tmp *= vec;
        psi += tmp;
        vec *= inv;
        vec += psi;
    }

    ValueType sum = vec.sum();

    std::cout << "Result is " << sum << std::endl;
}


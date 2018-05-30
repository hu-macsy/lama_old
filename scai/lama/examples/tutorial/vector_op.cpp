/**
 * @file vector.cpp
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
 * @brief vector.cpp
 * @author
 * @date 17.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;
using utilskernel::HArrayUtils;

typedef DefaultReal ValueType;

void methods2()
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10;

    auto x = fill<DenseVector<ValueType>>( n, 1, ctx );
    auto y = fill<DenseVector<ValueType>>( n, 2, ctx );

    x[0] = 0.5;
    y[1] = x[0] * 1.0 - 0.5 * y[0];

    x += 1.0;
    y -= 1.3;
    y *= 1.5;
    x /= 0.7;

    x += y;
    y -= x;
    x /= y;
    x *= y;

    y += x *= 2;

    // UnaryOp operations

    // x = 1 / x;            // x[i] = 1.0 / x[i]
    x.unaryOp( x, common::UnaryOp::RECIPROCAL );
    y = conj( y );        // y[i] = conj( y[i] )
    x = log( x );
    y = floor( y );
    x = ceil( x );
    x = sign( x );
    x = sqrt( x );
    x = sin( x );
    x = cos( x );
    x = tan( x );
    x = atan( x);
    x = pow( 2.0, x );
    // x.binaryOp( ValueType( 2 ), common::BinaryOp::POW, x );  // x[i] = 2.0 ** x[i] 
    x = pow( x, 2.0 );
    // y.binaryOp( y, common::BinaryOp::POW, ValueType( 2 ) );  // x[i] = x[i] ** 2.0
    x = pow( y, x );
    // x.binaryOp( y, common::BinaryOp::POW, x );  // x[i] = y[i] ** x[i]
    x = pow( x, y );
    // y.binaryOp( y, common::BinaryOp::POW, x );  // y[i] = y[i] ** x[i]

    ValueType s;

    s = x.sum();
    s = x.min();
    s = x.max();

    s = x.l1Norm();
    s = x.l2Norm();
    s = y.maxNorm();
   
    s = x.dotProduct( y );
    s = x.maxDiffNorm( y );

    std::cout << "max( abs( x - y ) ) = " << s << std::endl;

    x = 2 * x * y ;
    x = 2 * x - y;
}

int main()

{
    IndexType n = 10;

    auto xD = fill<DenseVector<ValueType>>( n, 1 );
    auto xS = fill<SparseVector<ValueType>>( n, 1 );

    std::cout << "DenseVector = " << xS << std::endl;

    xD = xS;

    ValueType s = xS[0];
    ValueType d = xD[0];

    std::cout << "xD[0] = " << d << ", xS[0] = " << s << std::endl;

    xD = 1.0;

    IndexType rawNonZeroIndexes1[] = { 0, 5, 6 };
    ValueType rawNonZeroValues1[] = { 0.5, 0.6, 0.7 };

    xS.setSparseValues( 3, rawNonZeroIndexes1, rawNonZeroValues1, ValueType( 1 ) );
 
    xD *= xS;  // be careful, sets most of xD to 0

    IndexType rawNonZeroIndexes2[] = { 0, 4, 6 };
    ValueType rawNonZeroValues2[] = { 0.5, 0.4, 0.9 };

    SparseVector<ValueType> xS1( n, 3, rawNonZeroIndexes1, rawNonZeroValues1, ValueType( 1 ) );
    SparseVector<ValueType> xS2( n, 3, rawNonZeroIndexes2, rawNonZeroValues2, ValueType( 1 ) );

    auto xD1 = convert<DenseVector<ValueType>>( xS1 );
    auto xD2 = convert<DenseVector<ValueType>>( xS2 );

    std::cout << "Max norm xS1 = " << xS1.maxNorm() << ", xD1 = " << xD1.maxNorm() << std::endl;
    std::cout << "L1 norm xS1 = " << xS1.l1Norm() << ", xD1 = " << xD1.l1Norm() << std::endl;
    std::cout << "L2 norm xS1 = " << xS1.l2Norm() << ", xD1 = " << xD1.l2Norm() << std::endl;
    std::cout << "Dotp xS1 * xS2 = " << xS1.dotProduct( xS2 ) << ", xD1 * xD2 = " << xD1.dotProduct( xD2 ) << std::endl;
}


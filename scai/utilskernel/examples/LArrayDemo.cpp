/**
 * @file LArrayDemo.cpp
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
 * @brief Example program that uses operations for LArrays.
 * @author
 * @date 17.03.2017
 */

#include <scai/utilskernel/LArray.hpp>

using namespace scai;
using utilskernel::LArray;

int main()
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10000;

    LArray<double> x( n, 1.0, ctx );
    LArray<double> y( n, 2.0, ctx );
 
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

    // unary operations

    x.invert();      // x[i] = 1.0 / x[i]
    y.conj();        // y[i] = conj( y[i] )
    x.log();
    y.floor();
    x.ceil();
    x.sqrt();
    x.sin();
    x.cos();
    x.tan();
    x.atan();
    x.powBase( 2.0 );  // x[i] = 2.0 ** x[i] 
    y.powExp( 2.0 );   // x[i] = x[i] ** 2.0
    x.powBase( y );    // x[i] = y[i] ** x[i]
    y.powExp( x );     // y[i] = y[i] ** x[i]

    double s;

    s = x.sum();
    s = x.min();
    s = x.max();

    s = x.l1Norm();
    s = x.l2Norm();
    s = y.maxNorm();
   
    s = x.dotProduct( y );
    s = x.maxDiffNorm( y );

    std::cout << "s = " << s << std::endl;
}

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
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>

#include <scai/utilskernel/LArray.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;
using utilskernel::LArray;

void methods1()
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10;

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
}

void methods2()
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10;

    DenseVector<double> x( n, 1.0, ctx );
    DenseVector<double> y( n, 2.0, ctx );

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

    Scalar s;

    s = x.sum();
    s = x.min();
    s = x.max();

    s = x.l1Norm();
    s = x.l2Norm();
    s = y.maxNorm();
   
    s = x.dotProduct( y );
    s = x.maxDiffNorm( y );

    x = 2 * x * y ;
    x = 2 * x - y;
}

int main()

{
    IndexType n = 10;

    DenseVector<double> xD( n, 1 );

    SparseVector<double> xS( n, 1 );

    std::cout << "DenseVector = " << xS << std::endl;

    xD = xS;

    Scalar s = xS[0];
    Scalar d = xD[0];

    std::cout << "xD[0] = " << d << ", xS[0] = " << s << std::endl;

    xD = 1.0;

    IndexType rawNonZeroIndexes1[] = { 0, 5, 6 };
    double rawNonZeroValues1[] = { 0.5, 0.6, 0.7 };

    xS.setSparseValues( 3, rawNonZeroIndexes1, rawNonZeroValues1, 1.0 );
 
    xD *= xS;  // be careful, sets most of xD to 0

    IndexType rawNonZeroIndexes2[] = { 0, 4, 6 };
    double rawNonZeroValues2[] = { 0.5, 0.4, 0.9 };

    SparseVector<double> xS1( n, 3, rawNonZeroIndexes1, rawNonZeroValues1, 1.0 );
    SparseVector<double> xS2( n, 3, rawNonZeroIndexes2, rawNonZeroValues2, 1.0 );

    DenseVector<double> xD1( xS1 );
    DenseVector<double> xD2( xS2 );

    std::cout << "Max norm xS1 = " << xS1.maxNorm() << ", xD1 = " << xD1.maxNorm() << std::endl;
    std::cout << "L1 norm xS1 = " << xS1.l1Norm() << ", xD1 = " << xD1.l1Norm() << std::endl;
    std::cout << "L2 norm xS1 = " << xS1.l2Norm() << ", xD1 = " << xD1.l2Norm() << std::endl;
    std::cout << "Dotp xS1 * xS2 = " << xS1.dotProduct( xS2 ) << ", xD1 * xD2 = " << xD1.dotProduct( xD2 ) << std::endl;
}


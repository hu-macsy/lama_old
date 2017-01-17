/**
 * @file lama/examples/tutorial/vector_exp.cpp
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
 * @brief ToDo: Missing description in ./lama/examples/tutorial/vector_exp.cpp
 * @author Lauretta Schubert
 * @date 14.08.2015
 */
#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/expression/all.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::lama;

int main()

{
    //
    // Define the ValueType used for the vector
    //
    typedef RealType ValueType;
    //
    // Vector expressions
    //
    DenseVector<ValueType> x, y, z;
    x = 1.0;
    y = 2.0;
    z = x + y;
    z = x * 2.0 + y;
    z = 2.0 * x + y;
    z = x + 1.0 * y;
    z = x + y * 1.0;
    z = y * 2.0;
    z = y / 2.0;
    z += x;
    z += 2.0 * x;
    z += x * 2.0;
    z -= x;
    z -= 2.0 * x;
    z -= x * 2.0;
    z *= 3.0;
    z /= 1.5;
    DenseVector<ValueType> tmp1 ( x + y );
    DenseVector<ValueType> tmp2 ( x * 2.0 + y );
    DenseVector<ValueType> tmp3 ( 2.0 * x + y );
    DenseVector<ValueType> tmp4 ( x + 1.0 * y );
    DenseVector<ValueType> tmp5 ( x + y * 1.0 );
    DenseVector<ValueType> tmp6 ( y * 2.0 );
    DenseVector<ValueType> tmp7 ( y / 2.0 );
    //
    // Matrix vector expressions
    //
    DenseMatrix<ValueType> A;
    z = A * x + 2.0 * y;
    z = A * x + y;
    z = A * x ;
    z = A * x + y * 2.0;
    DenseVector<ValueType> v1( A * x + 2.0 * y );
    DenseVector<ValueType> v2( A * x + y );
    DenseVector<ValueType> v3( A * x );
    DenseVector<ValueType> v4( A * x + y * 2.0 );
    z = 2.0 * A * x + 2.0 * y;
    z = 2.0 * A * x + y;
    z = 2.0 * A * x ;
    z = 2.0 * A * x + y * 2.0;
    DenseVector<ValueType> v5( 2.0 * A * x + 2.0 * y );
    DenseVector<ValueType> v6( 2.0 * A * x + y );
    DenseVector<ValueType> v7( 2.0 * A * x );
    DenseVector<ValueType> v8( 2.0 * A * x + y * 2.0 );
    z += A * x;
    z += 2.0 * A * x;
    z -= A * x;
    z -= 2.0 * A * x;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}


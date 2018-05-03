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

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;

int main()

{
    const IndexType N = 100;
    //
    // Define the ValueType used for the vector
    //
    typedef DefaultReal ValueType;
    //
    // Vector expressions
    //
    auto x = fill<DenseVector<ValueType>>( N, 1 );
    auto y = fill<DenseVector<ValueType>>( N, 2 );
    auto z = fill<DenseVector<ValueType>>( N, 3 );

    x = 1;
    y = 2;
    z = -x + y;
    z = - ( x * 2 + y );
    z = 2 * x + y;
    z = -y;
    z = 2 * y;
    z = 2 * ( 3 * y );
    z = 2 * (-y);
    z = x + 1 * (-y);
    z = x + 1 * y * 1;
    z = - ( y * 2 );
    z = - ( y / 2 );
    z += x;
    z += 2 * x;
    z += x * 2;
    z -= x;
    z -= 2 * x;
    z -= x * 2;
    z *= 3;
    z /= 1.5;

    auto tmp1 = eval<DenseVector<ValueType>>( x + y );
    auto tmp2 = eval<DenseVector<ValueType>>( x * 2 + y );
    auto tmp3 = eval<DenseVector<ValueType>>( 2 * x + y );
    auto tmp4 = eval<DenseVector<ValueType>>( x + 1 * y );
    auto tmp5 = eval<DenseVector<ValueType>>( x + y * 1 );
    auto tmp6 = eval<DenseVector<ValueType>>( y * 2 );
    auto tmp7 = eval<DenseVector<ValueType>>( y / 2 );
    //
    // _Matrix vector expressions
    //
    auto A = identity<DenseMatrix<ValueType>>( N );

    z = A * x + 2 * y;
    z = A * x + y;
    z = A * x ;
    z = A * x + y * 2;

    auto v1 = eval<DenseVector<ValueType>>( A * x + 2 * y );
    auto v2 = eval<DenseVector<ValueType>>( A * x + y );
    auto v3 = eval<DenseVector<ValueType>>( A * x );
    auto v4 = eval<DenseVector<ValueType>>( A * x + y * 2 );

    z = 2 * A * x + 2 * y;
    z = 2 * A * x + y;
    z = 2 * A * x ;
    z = 2 * A * x + y * 2;

    auto v5 = eval<DenseVector<ValueType>>( 2 * A * x + 2 * y );
    auto v6 = eval<DenseVector<ValueType>>( 2 * A * x + y );
    auto v7 = eval<DenseVector<ValueType>>( 2 * A * x );
    auto v8 = eval<DenseVector<ValueType>>( 2 * A * x + y * 2 );

    z += A * x;
    z += 2 * A * x;
    z -= A * x;
    z -= 2 * A * x;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}


/**
 * @file lama/examples/tutorial/vector_exp.cpp
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
    auto x = denseVectorFill<ValueType>( N, 1 );
    auto y = denseVectorFill<ValueType>( N, 2 );
    auto z = denseVectorFill<ValueType>( N, 3 );

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

    auto tmp1 = denseVectorEval( x + y );
    auto tmp2 = denseVectorEval( x * 2 + y );
    auto tmp3 = denseVectorEval( 2 * x + y );
    auto tmp4 = denseVectorEval( x + 1 * y );
    auto tmp5 = denseVectorEval( x + y * 1 );
    auto tmp6 = denseVectorEval( y * 2 );
    auto tmp7 = denseVectorEval( y / 2 );
    //
    // _Matrix vector expressions
    //
    auto A = identity<DenseMatrix<ValueType>>( N );

    z = A * x + 2 * y;
    z = A * x + y;
    z = A * x ;
    z = A * x + y * 2;

    auto v1 = denseVectorEval( A * x + 2 * y );
    auto v2 = denseVectorEval( A * x + y );
    auto v3 = denseVectorEval( A * x );
    auto v4 = denseVectorEval( A * x + y * 2 );

    z = 2 * A * x + 2 * y;
    z = 2 * A * x + y;
    z = 2 * A * x ;
    z = 2 * A * x + y * 2;

    auto v5 = denseVectorEval( 2 * A * x + 2 * y );
    auto v6 = denseVectorEval( 2 * A * x + y );
    auto v7 = denseVectorEval( 2 * A * x );
    auto v8 = denseVectorEval( 2 * A * x + y * 2 );

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


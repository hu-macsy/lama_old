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

#include <scai/common/Complex.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/expression/CastVectorExpression.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace lama;

int main()

{
    /** Take default real type for this example. */

    typedef float ValueType;

    typedef common::Complex<float> ComplexType;

    DenseVector<ValueType> x = linearValuesVector( 10, ValueType( 1 ), ValueType( 0.5 ) );
    DenseVector<ValueType> y = linearValuesVector( 10, ValueType( 20 ), ValueType( -0.5 ) );

    DenseVector<ComplexType> z;
    z = cmplx( x, y );

    z.writeToFile( "z.txt" );

    x = real( z );
    y = imag( z );

    x.writeToFile( "x.txt" );
    y.writeToFile( "y.txt" );

    return EXIT_SUCCESS;
}


/**
 * @file vector.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief vector.cpp
 * @author
 * @date 17.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::lama;

int main()

{
    /** Take default real type for this example. */

    typedef RealType ValueType;

    Scalar singleValue( 2.0 );

    //
    // Create a DenseVector out of a simple c array
    //
    const ValueType inputData[] = { 1.0, 2.0, 3.0, 4.0 };

    DenseVector<ValueType> sequenceOfValues( 4, inputData );

    //
    // scale vector
    //
    sequenceOfValues = singleValue * sequenceOfValues;

    //
    // print vector to file vector.frm/.vec (SAMG format)
    //
    sequenceOfValues.writeToFile( "vector", File::SAMG_FORMAT );

    std::cout << "DenseVector is written to 'vector.frm/.vec'" << std::endl;

    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

    return EXIT_SUCCESS;
}


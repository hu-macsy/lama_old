/**
 * @file vector.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief vector.cpp
 * @author
 * @date 17.05.2013
 * @since 1.0.0
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
    //
    // Define the ValueType used for the vector
    //
    typedef double ValueType;

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
    sequenceOfValues.writeToFile( "vector", File::FORMATTED );

    std::cout << "DenseVector is written to 'vector.frm/.vec'" << std::endl;

    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

    return EXIT_SUCCESS;
}


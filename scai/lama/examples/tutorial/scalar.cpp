/**
 * @file lama/examples/tutorial/scalar.cpp
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
 * @brief ToDo: Missing description in ./lama/examples/tutorial/scalar.cpp
 * @author Lauretta Schubert
 * @date 14.08.2015
 */
#include <scai/lama.hpp>

#include <scai/lama/Scalar.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::lama;

int main()

{
    //
    // Create scalars
    //
    Scalar a( 1 );    // a scalar of index type
    Scalar b( 2.5f ); // a scalar of type float
    Scalar c( 0.0 );  // a scalar representing zero
    //
    // binary operators
    //
    c = a + b;
    c = a - b;
    c = a * b;
    c = a / b;
    a += b;
    a -= b;
    a *= b;
    a /= b;
    //
    // UnaryOp operator '-'
    //
    c = -c;
    //
    //relational operators
    //
    bool boolean = ( a == b );
    boolean = ( a != b );
    boolean = ( a < b );
    boolean = ( a > b );
//    boolean = ( a <= b );
//    boolean = ( a >= b );
//    std::cout << "a >= B : " << boolean << std::endl;
    std::cout << "a > B : " << boolean << std::endl;
    //
    // math functions
    //
    c = max( a, b );
    c = min( a, b );
    c = abs( a );
    c = sqrt( a );
    //
    // output operator
    //
    std::cout << "my Scalar is: " << a << std::endl;
    // alternative: getValue
//    printf("my Scalar is: %d\n", a.getValue<int>() );
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}

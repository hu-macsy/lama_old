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
    // unary operator '-'
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

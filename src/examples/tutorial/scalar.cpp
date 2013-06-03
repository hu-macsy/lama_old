#include <lama.hpp>

#include <lama/Scalar.hpp>

#include <iostream>

using namespace lama;

int main()

{
    Scalar a(1);    // a scalar of index type
    Scalar b(2.5f); // a scalar of type float
    Scalar c(0.0);  // a scalar representing zero

    // binary operators
    c = a + b;
    c = a - b;
    c = a * b;
    c = a / b;

    a += b;
    a -= b;
    a *= b;
    a /= b;

    // unary operator '-'
    c = -c;

    //relational operators
    bool boolean = ( a == b );
    boolean = ( a != b );
    boolean = ( a < b );
    boolean = ( a > b );
    boolean = ( a <= b );
    boolean = ( a >= b );

    // math functions
    c = max( a, b );
    c = min( a, b );
    c = abs( a );
    c = sqrt( a );

    // output operator
    std::cout << "my Scalar is: " << a << std::endl;

    // getValue
    printf("my Scalar is: %d\n", a.getValue<int>() );

    return 0;
}

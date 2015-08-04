#include <lama.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>

#include <iostream>
#include <stdlib.h>

using namespace lama;

int main()

{
    //
    // Define the ValueType used for the vector
    //
    typedef double ValueType;

    //
    // Create a DenseVector of size 8 with value 1.1 in each row
    //
    IndexType size = 8;
    DenseVector<ValueType> v( size, 1.1 );

    //
    // Getting the L1 norm of the vector and printing out the scalar
    //
    Scalar s = v.l1Norm();
    std::cout << "L1 norm of v = " << s.getValue<ValueType>() << std::endl;

    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

    return EXIT_SUCCESS;
}

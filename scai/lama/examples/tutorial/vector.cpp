#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>

#include <iostream>
#include <stdlib.h>

using namespace lama;

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


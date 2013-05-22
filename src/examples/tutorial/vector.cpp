#include <lama.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/expression/all.hpp>

#include <iostream>

using namespace lama;

int main()

{
    typedef double ValueType;

    Scalar singleValue( 2.0 );

    const ValueType inputData[] = { 1.0, 2.0, 3.0, 4.0 };

    DenseVector<ValueType> sequenceOfValues( 4, inputData );

    sequenceOfValues = singleValue * sequenceOfValues;

    sequenceOfValues.writeToFile( "vector.txt" , File::FORMATTED );
}


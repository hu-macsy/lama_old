#include <lama.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/matrix/DenseMatrix.hpp>
#include <lama/expression/all.hpp>

#include <iostream>

using namespace lama;

int main()

{
    typedef double ValueType;

    DenseVector<ValueType> x, y, z;

    x = 1.0;
    y = 2.0;

    z = x + y;
    z = x * 2.0 + y;
    z = 2.0 * x + y;
    z = x + 1.0 * y;
    z = x + y * 1.0;
    z = y * 2.0;
    z = y / 2.0;
    z += x;
    z += 2.0 * x;
    z += x * 2.0;
    z -= x;
    z -= 2.0 * x;
    z -= x * 2.0;
    z *= 3.0;
    z /= 1.5;

    DenseVector<ValueType> tmp1 ( x + y );
    DenseVector<ValueType> tmp2 ( x * 2.0 + y );
    DenseVector<ValueType> tmp3 ( 2.0 * x + y );
    DenseVector<ValueType> tmp4 ( x + 1.0 * y );
    DenseVector<ValueType> tmp5 ( x + y * 1.0 );
    DenseVector<ValueType> tmp6 ( y * 2.0 );
    DenseVector<ValueType> tmp7 ( y / 2.0 );

    DenseMatrix<ValueType> A;

    z = A * x + 2.0 * y;
    z = A * x + y;
    z = A * x ;
    z = A * x + y * 2.0;

    DenseVector<ValueType> v1( A * x + 2.0 * y );
    DenseVector<ValueType> v2( A * x + y );
    DenseVector<ValueType> v3( A * x );
    DenseVector<ValueType> v4( A * x + y * 2.0 );

    z = 2.0 * A * x + 2.0 * y;
    z = 2.0 * A * x + y;
    z = 2.0 * A * x ;
    z = 2.0 * A * x + y * 2.0;

    DenseVector<ValueType> v5( 2.0 * A * x + 2.0 * y );
    DenseVector<ValueType> v6( 2.0 * A * x + y );
    DenseVector<ValueType> v7( 2.0 * A * x );
    DenseVector<ValueType> v8( 2.0 * A * x + y * 2.0 );

    z += A * x;
    z += 2.0 * A * x;
    z -= A * x;
    z -= 2.0 * A * x;
}


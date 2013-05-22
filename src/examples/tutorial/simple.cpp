#include <lama.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>

#include <iostream>

using namespace lama;

int main()

{
    IndexType size = 8;
    DenseVector<double> v( size, 1.1 );
    Scalar s = v.l1Norm();
    std::cout << "L1 norm of v = " << s.getValue<double>() << std::endl;
    return 0;
}


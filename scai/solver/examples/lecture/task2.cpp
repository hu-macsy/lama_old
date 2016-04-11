
//Solution of task 2:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/all.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/tracing.hpp>

#include <iostream>

using namespace scai::lama;
using namespace scai::hmemo;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    int maxIter = 10;   // maximal number of iterations

    if ( argc > 2 ) 
    {
       sscanf( argv[2], "%d", &maxIter );
    }

    CSRSparseMatrix<double> A( argv[1] );
    std::cout << "Read matrix A : " << A << std::endl;
    IndexType size = A.getNumRows();

    DenseVector<double> b( size, 0 );

    {
        WriteAccess<double> writeB( b.getLocalValues() );

        for ( IndexType i = 0; i < size; ++i )
        {
            // writeB[i] = 1;
            writeB[i] = double( i + 1 );
        }
    }

    std::cout << "Vector b : " << b << std::endl;
    DenseVector<double> x( size , 0.0 );
    std::cout << "Vector x : " << x << std::endl;

    // d = r = b - A * x
    // help = A * x;

    DenseVector<double> r = b - A * x;
    DenseVector<double> d = r;
    Scalar rOld = r.dotProduct( r );
    Scalar eps = 0.00001;

    L2Norm norm;

    for ( int k = 0 ; k < maxIter and norm(r) > eps; k++ )
    {
        DenseVector<double> z = A * d;
        Scalar alpha = rOld / d.dotProduct( z );
        x = x + alpha * d;
        r = r - alpha * z;
        Scalar rNew = r.dotProduct( r );
        Scalar beta = rNew / rOld;
        d = r + beta * d;
        rOld = rNew;

        Scalar rnorm = norm( r );
        std::cout << "Iter k = " << k << " : norm( r ) = " << rnorm.getValue<double>() << std::endl;
    }

    return 0;
}


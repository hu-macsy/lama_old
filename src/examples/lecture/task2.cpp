
//Solution of task 2:

#include <lama.hpp>

#include <lama/storage/SparseAssemblyStorage.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/DenseVector.hpp>

#include <lama/expression/all.hpp>
#include <lama/norm/L2Norm.hpp>

#include <tracing/tracing.hpp>

#include <iostream>

using namespace lama;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    CSRSparseMatrix<double> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();
    DenseVector<double> rhs( size , 0.0 );
    HostWriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    DenseVector<double> solution( size , 0.0 );
    std::cout << "Vector solution : " << solution << std::endl;
    Scalar rNew;
    Scalar rOld;
    Scalar alpha;
    Scalar beta;
    Scalar eps = 0.00001;
    DenseVector<double> r ( size, 0.0 );
    DenseVector<double> help ( size, 0.0 );
    DenseVector<double> d ( size, 0.0 );
    DenseVector<double> Ad( size, 0.0 );
    // d = r = solution - m*rhs
    help = m * solution;
    r = rhs - help;
    d = r;
    rOld = r.dotProduct( r );
    L2Norm norm;

    for ( int k = 1 ; k < 10 ; k++ )
    {
        Scalar rnorm = norm( r );

        std::cout << "Iter k = " << k << " : norm( r ) = " << rnorm.getValue<double>() << std::endl;

        if ( norm( r ) < eps )
        {
            break;
        }
        else
        {
            Ad = m * d;
            alpha = rOld / ( d.dotProduct( Ad ) );
            solution = solution + ( alpha * d ) ;
            r = r - alpha * Ad;
            rNew = r.dotProduct( r );
            beta = rNew / rOld;
            d = r + beta * d;
            rOld = rNew;
        }
    }

    return 0;
}


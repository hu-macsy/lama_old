
//Solution of task 2a:

#include <lama.hpp>

#include <lama/storage/SparseAssemblyStorage.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/DenseVector.hpp>
#include <tracing.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/norm/L2Norm.hpp>

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
    WriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    DenseVector<double> solution( size , 0.0 );
    std::cout << "Vector solution : " << solution << std::endl;
    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    CG cgSolver( "CGTestSolver" );
    cgSolver.setStoppingCriterion( rt );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );
    std::cout << "The solution is: ";
    ReadAccess<double> hra( solution.getLocalValues() );

    for ( int i = 0; i < solution.size(); i++ )
    {
        std::cout << hra[i] << " ";
    }

    hra.release();
    return 0;
}


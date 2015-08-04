
//Solution of task 4:

#include <lama.hpp>

#include <lama/storage/SparseAssemblyStorage.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/DenseVector.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/criteria/IterationCount.hpp>

#include <lama/distribution/BlockDistribution.hpp>

#include <tracing.hpp>

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

    CommunicatorPtr comm( Communicator::get( "MPI" ) );
    DistributionPtr dist( new BlockDistribution( size, comm ) );
    m.redistribute( dist, dist );

    DenseVector<double> rhs( size , 0.0 );
    WriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    rhs.redistribute( dist );
    DenseVector<double> solution( size, 0.0 );
    solution.redistribute( dist );
    std::cout << "Vector solution : " << solution << std::endl;
    CG cgSolver( "CGTestSolver" );
    CriterionPtr criterion( new IterationCount( 10 ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );

    std::cout << "The solution is: ";
    for ( int i = 0; i < solution.size(); ++i )
    {
        std::cout << solution.getValue( i ) << " ";
    }
    std::cout << std::endl;
    return 0;
}


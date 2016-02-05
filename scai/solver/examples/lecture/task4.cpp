
//Solution of task 4:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/distribution/BlockDistribution.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/tracing.hpp>

#include <iostream>

using namespace scai;
using namespace lama;
using namespace solver;
using namespace hmemo;

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

    CommunicatorPtr comm( Communicator::getCommunicator( scai::lama::communicator::MPI ) );
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


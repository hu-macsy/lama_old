
//Solution of task 5:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Communicator.hpp>

#include <scai/lama/solver/CG.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>

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

    CSRSparseMatrix<double> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();

    ContextPtr cudaContext = Context::getContextPtr( scai::hmemo::context::CUDA, 0 ); 
    m.setContext( cudaContext );

    DenseVector<double> rhs( size , 0.0 );
    WriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    rhs.setContext( cudaContext );
    DenseVector<double> solution( size, 0.0 );
    solution.setContext( cudaContext );
    std::cout << "Vector solution : " << solution << std::endl;
    CG cgSolver( "CGTestSolver" );
    CriterionPtr criterion( new IterationCount ( 10 ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );

    return 0;
}



//Solution of task 5:

#include <lama.hpp>

#include <lama/storage/SparseAssemblyStorage.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/criteria/IterationCount.hpp>

#include <lama/distribution/BlockDistribution.hpp>

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

    lama::ContextPtr cudaContext = ContextFactory::getContext( Context::CUDA, 0 ); 
    m.setContext( cudaContext );

    DenseVector<double> rhs( size , 0.0 );
    HostWriteAccess<double> hwarhs( rhs.getLocalValues() );

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


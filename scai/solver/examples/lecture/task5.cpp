
//Solution of task 5:

#include <scai/lama.hpp>

#include <scai/lama/storage/SparseAssemblyStorage.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Communicator.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <iostream>

using namespace scai;
using namespace scai::lama;
using namespace scai::solver;
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

    ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA, 0 ); 
    m.setContextPtr( cudaContext );

    DenseVector<double> rhs( size , 0.0 );
    WriteAccess<double> hwarhs( rhs.getLocalValues() );

    for ( IndexType i = 0; i < size; ++i )
    {
        hwarhs[i] = double( i + 1 );
    }

    std::cout << "Vector rhs : " << rhs << std::endl;
    hwarhs.release();
    rhs.setContextPtr( cudaContext );
    DenseVector<double> solution( size, 0.0 );
    solution.setContextPtr( cudaContext );
    std::cout << "Vector solution : " << solution << std::endl;
    CG cgSolver( "CGTestSolver" );
    CriterionPtr criterion( new IterationCount ( 10 ) );
    cgSolver.setStoppingCriterion( criterion );
    cgSolver.initialize( m );
    cgSolver.solve( solution, rhs );

    return 0;
}


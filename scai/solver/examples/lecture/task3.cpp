
//Solution of task 3:

#include <scai/lama.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>

#include <iostream>

using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace solver;

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "No input file specified" << std::endl;
        exit( -1 );
    }

    //Read a sparse matrix from the passed input file
    CSRSparseMatrix<double> m( argv[1] );
    std::cout << "Read matrix m : " << m << std::endl;
    IndexType size = m.getNumRows();

    //Create rhs vector
    DenseVector<double> rhs( size , 0.0 );
    std::cout << "Vector rhs : " << rhs << std::endl;

    //Create solution vector
    DenseVector<double> solution( size , 1.0 );
    std::cout << "Vector solution : " << solution << std::endl;

    //Compute the rhs that fits our solution to be able to calculate the error later
    rhs = m * solution;

    //Reset solution to zero so that there is something to solve
    solution = 0.0;

    common::shared_ptr<Timer> timer( new Timer() );
    LoggerPtr logger( new CommonLogger("CGLogger", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly, timer) );

    //Create a CG solver using the Constructor with logger as an argument
    CG cgSolver( "CGTestSolver", logger );

    //Create a stopping criterion for the iterative solver cgSolver
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr criterion( new ResidualThreshold( norm, 1E-8, ResidualThreshold::Absolute ) );
    cgSolver.setStoppingCriterion( criterion );

    //Initialize the solver
    cgSolver.initialize( m );
    //Solve m * solution = rhs
    cgSolver.solve( solution, rhs );

    //calculate the error and its L2-Norm
    DenseVector<double> error( size, 1.0 );

    error = error - solution;

    std::cout<<"L2-Norm of error is "<<l2Norm( error )<<std::endl;

    return 0;
}


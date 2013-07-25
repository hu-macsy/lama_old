#include "Config.hpp"

#include <lama/solver/BiCG.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/DenseVector.hpp>

#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace lama;

/** ValueType is the type used for matrix and vector elements. */

typedef float ValueType;

int main( int argc, char* argv[] )
{
    // init

    Config config;

    const char* matrixFilename;
    const char* rhsFilename;
    const char* solutionFilename;
    
//    if ( argc >= 4 )
//    {
        matrixFilename = argv[1];
//        rhsFilename = argv[2];
//        solutionFilename = argv[3];
//
//        for ( int i = 4; i < argc; ++i )
//        {
//            config.setArg( argv[i] );
//        }
//    }
//    else
//    {
//        cout << "Usage: " << argv[0] << " <matrix.frm> <rhs.frv> <solution.frv>" << endl;
//        exit( 1 );
//    }

    CSRSparseMatrix<ValueType> csrMatrix( matrixFilename );
    DenseVector<ValueType> rhs( 675, 1.0);//rhsFilename );
    DenseVector<ValueType> solution( 675, 0.0); //solutionFilename );

    // BiCG

    LoggerPtr slogger(
        new CommonLogger( "<BiCG>: ", LogLevel::completeInformation, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new Timer() ) ) );

    BiCG BiCGSolver( "BiCGTestSolver", slogger );

    IndexType expectedIterations = 50;
    CriterionPtr criterion( new IterationCount( expectedIterations ) );

    BiCGSolver.setStoppingCriterion( criterion );
    BiCGSolver.initialize( csrMatrix );

    BiCGSolver.solve( solution, rhs );

}

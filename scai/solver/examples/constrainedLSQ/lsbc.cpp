/**
 * @file lsbc.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Example of least square problem with boundary conditions
 * @author Thomas Brandes, Andreas Borgen Longva
 * @date 21.07.2017
 */

#include "ConstrainedLeastSquares.hpp"
#include "MatrixWithT.hpp"

#include <scai/tracing.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/partitioning/Partitioning.hpp>

#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;
using namespace solver;
using namespace partitioning;

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.driver" )
    
    common::Settings::parseArgs( argc, argv );

    if ( argc < 5 )
    {
        std::cout << "Please call " << argv[0] 
                  << " <filename_A> <filename_b> <filename_lb> <filename_ub> [ <filename_x> ]" << std::endl;
        return -1;
    }

    CSRSparseMatrix<double> A( argv[1] );
    DenseVector<double> b ( argv[2] );
    DenseVector<double> lb ( argv[3] );
    DenseVector<double> ub ( argv[4] );

    std::cout << "A = " << A << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "lb = " << lb << std::endl;
    std::cout << "ub = " << ub << std::endl;

    if ( false )
    {
        A.writeToFile( "A.mtx" );
        b.writeToFile( "b.mtx" );
        lb.writeToFile( "lb.mtx" );
        ub.writeToFile( "ub.mtx" );
    }

    // take context as specified by SCAI_CONTEXT

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    A.setContextPtr( ctx );
    A.setCommunicationKind( _Matrix::SYNCHRONOUS );
    b.setContextPtr( ctx );
    ub.setContextPtr( ctx );
    lb.setContextPtr( ctx );

    std::string partitioningKind = "BLOCK";

    common::Settings::getEnvironment( partitioningKind, "SCAI_PARTITIONING" );

    if ( !Partitioning::canCreate( partitioningKind ) )
    {
        std::cout << "ERROR: Partitioning kind = " << partitioningKind << " not supported" << std::endl;

        std::vector<std::string> values;  // string is create type for the factory

        Partitioning::getCreateValues( values );

        std::cout << "Supported partitiong types are:";

        for ( size_t i = 0; i < values.size(); ++i )
        {
            std::cout << " " << values[i];
        }

        std::cout << std::endl;

        return -1;
    }

    PartitioningPtr thePartitioning( Partitioning::create( partitioningKind ) );
  
    thePartitioning->rectangularRedistribute( A, ( float ) 1.0 );

    dmemo::DistributionPtr rowDist = A.getRowDistributionPtr();
    dmemo::DistributionPtr colDist = A.getColDistributionPtr();

    lb.redistribute( colDist );
    ub.redistribute( colDist );
    b.redistribute( rowDist );

    DenseVector<double> x( ctx );

    _MatrixWithT Aopt( A );   // Allocate also a transposed matrix to optimize A' * x operations

    ConstrainedLeastSquares lsq( Aopt );

    // lsq.useTranspose();       // will matrixTimesVector instead ov vectorTimesMatrix

    lsq.setTolerance( 0.01 );
    lsq.setMaxIter( 50 );

    try 
    {
        lsq.solve( x, b, lb, ub );
    }
    catch ( common::Exception& ex )
    {
        std::cout << "Caught exception: " << ex.what() << std::endl;
        std::cout << "Stop execution." << std::endl;
        return 1;
    }

    DenseVector<double> residual( A * x - b );

    std::cout << "res norm = " << residual.l2Norm() << std::endl;

    if ( argc > 5 )
    {
        x.writeToFile( argv[5] );
        std::cout << "written solution to file " << argv[5] << std::endl;
    }
}

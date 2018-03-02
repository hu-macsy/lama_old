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
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#include <ipbcls/ConstrainedLeastSquares.hpp>

#include <scai/tracing.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/lama/matrix/MatrixWithT.hpp>

#include <scai/common/Settings.hpp>
#include <scai/partitioning/Partitioning.hpp>

#include <cstdlib>

using namespace scai;
using namespace lama;

/** Set row and column distribution of a matrix 
 *
 *  @param[in,out] determine a good distribution of sparse matrix A and redistribute it
 *
 *  Default: BLOCK distribution, but graph partitioning via METIS/SCOTCH, .. might be used
 *
 *  uses environment variables:
 *
 *    SCAI_PARTITIONING ( e.g. BLOCK, METIS, ... )
 *    SCAI_WEIGHT=w1,w2,...,wp    comma separated list for individual weights of processors
 */
void setDistribution( _Matrix& A )
{
    using namespace partitioning;
    using namespace common;

    std::string partitioningKind = "BLOCK";   // default setting

    std::string val;

    if ( Settings::getEnvironment( val, "SCAI_PARTITIONING" ) )
    {
        if ( Partitioning::canCreate( val ) )
        {
            partitioningKind = val;
        }
        else
        {
            std::cout << "ERROR: Partitioning kind = " << val << " not supported" << std::endl;
        }
    }

    Settings::setRank( A.getRowDistribution().getCommunicator().getRank() );

    float weight = 1.0f;   // default weight on this processor

    if ( Settings::getEnvironment( val, "SCAI_WEIGHT" ) )
    {
        weight = strtof( val.c_str(), NULL );
    }

    PartitioningPtr thePartitioning( Partitioning::create( partitioningKind ) );

    thePartitioning->rectangularRedistribute( A, weight );
}

int main( int argc, const char* argv[] )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    SCAI_REGION( "Main.driver" )
    
    common::Settings::parseArgs( argc, argv );

    if ( argc < 6 )
    {
        if ( comm->getRank() == 0 )
        {
            std::cout << "Please call " << argv[0] 
                      << " <tolerance> <filename_A> <filename_b> <filename_lb> <filename_ub> [ <filename_x> ]" << std::endl
                      << " where <tolerance> is a floating point number between 0 and 1. "<< std::endl;
        }

        return -1;
    }
    
    double tol = std::strtod(argv[1], NULL);
    if (tol == 0.0)
    {
        std::cout << "Invalid tolerance supplied." << std::endl;
        return -1;
    }

    auto A  = read<CSRSparseMatrix<double>>( argv[2] );
    auto b  = read<DenseVector<double>>( argv[3] );
    auto lb = read<DenseVector<double>>( argv[4] );
    auto ub = read<DenseVector<double>>( argv[5] );

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
    A.setCommunicationKind( SyncKind::SYNCHRONOUS );
    b.setContextPtr( ctx );
    ub.setContextPtr( ctx );
    lb.setContextPtr( ctx );

    setDistribution( A );   // redistribute A by using a partitioning class

    dmemo::DistributionPtr colDist( A.getColDistributionPtr() );
    dmemo::DistributionPtr rowDist( A.getRowDistributionPtr() );

    lb.redistribute( colDist );
    ub.redistribute( colDist );
    b.redistribute( rowDist );

    std::cout << "A = " << A << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "lb = " << lb << std::endl;
    std::cout << "ub = " << ub << std::endl;

    DenseVector<double> x( colDist, 0, ctx );

    MatrixWithT<double> Aopt( A );   // Allocate also a transposed matrix to optimize A' * x operations

    ConstrainedLeastSquares<double> lsq( Aopt );

    // lsq.setInnerSolverType( InnerSolverType::StandardCG );
    lsq.setInnerSolverType( InnerSolverType::NewtonStepCG );

    // lsq.useTranspose();       // will matrixTimesVector instead ov vectorTimesMatrix

    lsq.setObjectiveTolerance( tol );
    lsq.setMaxIter( 500 );

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

    auto residual = eval<DenseVector<double>>( A * x - b );

    std::cout << "res norm = " << residual.l2Norm() << std::endl;

    if ( argc > 6 )
    {
        x.writeToFile( argv[6] );
        std::cout << "written solution to file " << argv[6] << std::endl;
    }
}

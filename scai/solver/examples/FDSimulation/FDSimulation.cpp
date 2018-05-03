/**
 * @file solver/examples/FDSimulation/FDSimulation.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./solver/examples/FDSimulation/FDSimulation.cpp
 * @author Lauretta Schubert
 * @date 27.09.2016
 */
#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>

#include <iostream>
#include <memory>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration.hpp"

using namespace scai;

#define MASTER 0

#define HOST_PRINT( comm, msg )             \
    {                                           \
        PartitionId myRank = comm->getRank();   \
        if ( myRank == MASTER )                 \
        {                                       \
            std::cout << msg;                   \
        }                                       \
    }

/*
 *  routine for initializing derivative matrices (needed for velocity Updates vz, vx, vy)
 *  create csr data structure (sub)matrices (and replicate along matrix diagonal)
 */
template<typename ValueType>
void derivatives( lama::SparseMatrix<ValueType>& A,
                  lama::SparseMatrix<ValueType>& B,
                  lama::SparseMatrix<ValueType>& C,
                  IndexType NX,
                  IndexType NY,
                  IndexType NZ,
                  dmemo::DistributionPtr dist,
                  hmemo::ContextPtr ctx,
                  dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "derivatives" )

    // _Matrix A,B are created in 2 steps:
    //   1: create MatrixStorage in CSR format
    //   2: duplicate MatrixStorage along Diagonal till full _Matrix size is reached

    // for creating CSR MatrixStorage
    std::unique_ptr<lama::MatrixStorage<ValueType> > storageHelp( new lama::CSRStorage<ValueType>() );
    IndexType numValues;
    std::vector<IndexType> csrIA;
    std::vector<IndexType> csrJA;
    std::vector<ValueType> csrValues;

    // create matrix A
    numValues = NZ + ( NZ - 1 ); // diagonal element (NZ) + secondary diagonal elements (NZ - 1)
    csrIA.reserve( NZ + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );

    IndexType count = 0;
    IndexType size = NZ;
    csrIA.push_back( 0 );

    for ( IndexType i = 0; i < size; ++i )
    {
        for ( IndexType j = 0; j < size; ++j )
        {
            if ( i == j )
            {
                ++count;
                csrJA.push_back( j );
                csrValues.push_back( -1.0 );
            }

            if ( j - 1 == i )
            {
                ++count;
                csrJA.push_back( j );
                csrValues.push_back( 1.0 );
            }
        }

        csrIA.push_back( count );
    }

    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    lama::MatrixCreator::buildReplicatedDiag( A, *storageHelp, NX * NY ); // not distributed, need to redistribute afterwards
    A.redistribute( dist, dist );
    A.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix A finished\n" );

    csrIA.clear();
    csrJA.clear();
    csrValues.clear();

    // create matrix B
    numValues = NZ * NX + ( NZ * NX - ( NZ + 1 ) + 1 ); // diagonal element (NZ*NX) + secondary diagonal elements (NZ*NX - (NZ+1) + 1)
    csrIA.reserve( NZ + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );

    count = 0;
    size = NZ * NX;
    csrIA.push_back( 0 );

    for ( IndexType i = 0; i < size; ++i )
    {
        for ( IndexType j = 0; j < size; ++j )
        {
            if ( i == j )
            {
                ++count;
                csrJA.push_back( j );
                csrValues.push_back( -1.0 );
            }

            if ( j - NZ == i )
            {
                ++count;
                csrJA.push_back( j );
                csrValues.push_back( 1.0 );
            }
        }

        csrIA.push_back( count );
    }

    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    lama::MatrixCreator::buildReplicatedDiag( B, *storageHelp, NY ); // not distributed, need to redistribute afterwards
    B.redistribute( dist, dist );
    B.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix B finished\n" );

    csrIA.clear();
    csrJA.clear();
    csrValues.clear();

    // create matrix C
    // initialize by diagonals

    PartitionId myRank   = comm->getRank();
    PartitionId numRanks = comm->getSize();

    IndexType globalSize = dist->getGlobalSize();
    IndexType numDiagonals = 2;
    IndexType secondaryIndex = NZ * NX; // = remaining part of secondary diagonal
    IndexType numSecondary = globalSize - secondaryIndex;

    IndexType lb;
    IndexType ub;
    dmemo::BlockDistribution::getLocalRange( lb, ub, dist->getGlobalSize(), myRank, numRanks );

    size = dist->getLocalSize(); //getGlobalSize();

    IndexType myStart = lb;
    IndexType myEnd   = std::min( ub, numSecondary );
    IndexType mySize  = std::max( myEnd - myStart, IndexType( 0 ) );
    IndexType myRemaining = size - mySize;

    std::vector<IndexType> offsets;
    // distributed offset start with their lb
    offsets.push_back( lb + 0 );              // index main diagonal (starting with '0')
    offsets.push_back( lb + secondaryIndex ); // index secondary diagonal

    std::vector<ValueType> diagonals;
    diagonals.reserve( numDiagonals * size );
    // insert backwards so we can insert from the beginning of the vector
    diagonals.insert( diagonals.begin(), myRemaining, 0.0 ); // add ghost secondary diagonal (out of column bound)
    diagonals.insert( diagonals.begin(), mySize,      1.0 ); // add secondary diagonal
    diagonals.insert( diagonals.begin(), size,       -1.0 ); // insert main diagonal before secondary diagonal

    // Define a 'global' DIAStorage, does not copy the values 

    lama::DIAStorage<ValueType> diaStorage( size, globalSize,
                                            hmemo::HArrayRef<IndexType>( offsets ), 
                                            hmemo::HArrayRef<ValueType>( diagonals ) );

    C.assignLocal( diaStorage, dist );
    C.redistribute( dist, dist );      // compute halos

    C.setContextPtr( ctx );

    HOST_PRINT( comm, "Matrix C finished\n" );
}

/*
 *  routine for initializing all system matrices
 *  uses derivatives for initializing A, B, C
 */
template<typename ValueType>
void initializeMatrices( lama::SparseMatrix<ValueType>& A, lama::SparseMatrix<ValueType>& B, lama::SparseMatrix<ValueType>& C,
                         lama::Matrix<ValueType>& D, lama::Matrix<ValueType>& E, lama::Matrix<ValueType>& F, 
                         dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,
                         IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "initializeMatrices" )

    HOST_PRINT( comm, "Initialization of the matrices A,B,C,D,E,F...\n" );

    derivatives( A, B, C, NX, NY, NZ, dist, ctx, comm );

    D.setContextPtr( ctx );
    E.setContextPtr( ctx );
    F.setContextPtr( ctx );

    D.assignTranspose( A );
    D.scale( -1 );
    HOST_PRINT( comm, "Matrix D finished\n" );

    E.assignTranspose( B );
    E *= -1;
    HOST_PRINT( comm, "Matrix E finished\n" );

    F.assignTranspose( C );
    F *= -1;
    HOST_PRINT( comm, "Matrix F finished\n" );

    HOST_PRINT( comm, "Finished with initialization of the matrices!\n" );
}

/*
 *  routine for calculation the input source (Ricker wavelet)
 *
 *  MATLAB:
 *  t=0:DT:(NT*DT-DT);
 *  tau=pi*FC*(t-1.5/FC);
 *  signal=AMP*(1-2*tau.^2).*exp(-tau.^2);
 */
template<typename ValueType>
void sourceFunction( lama::DenseVector<ValueType>& source, IndexType FC, IndexType AMP, dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "sourceFunction" )

    HOST_PRINT( comm, "Calculate source signal...\n" );

    // this is for tau[i] = pi * FC * ( source[i] - 1.5/FC );

    auto help = lama::fill<lama::DenseVector<ValueType>>( source.size(), 1.5 / FC );
    auto tau  = lama::eval<lama::DenseVector<ValueType>>( source - help );

    tau *= M_PI * FC;

    // this is for source[i] = AMP * ( 1.0 - 2.0 * tau[i] * tau[i] * exp( -tau[i] * tau[i] ) );

    auto one = lama::fill<lama::DenseVector<ValueType>>( source.size(), 1 );
    help = tau * tau;
    tau = -help;
    tau = exp( tau );
    help = one - 2.0 * help;
    source = ValueType( AMP ) * help * tau;
}

/*
 *  routine doing NT time steps updating vX, vY, vZ, p
 *  with incoming source
 *  storing seismogram data
 */
template <typename ValueType>
void timesteps( lama::DenseVector<ValueType>& seismogram, lama::DenseVector<ValueType>& source, lama::DenseVector<ValueType>& p,
                lama::Vector<ValueType>& vX, lama::Vector<ValueType>& vY, lama::Vector<ValueType>& vZ,
                lama::Matrix<ValueType>& A, lama::Matrix<ValueType>& B, lama::Matrix<ValueType>& C, 
                lama::Matrix<ValueType>& D, lama::Matrix<ValueType>& E, lama::Matrix<ValueType>& F,
                ValueType v_factor, ValueType p_factor,
                IndexType NT, ValueType DH_INV, IndexType source_index, IndexType seismogram_index,
                dmemo::CommunicatorPtr comm, dmemo::DistributionPtr /*dist*/ )
{
    SCAI_REGION( "timestep" )

    for ( IndexType t = 0; t < NT; t++ )
    {
        if ( t % 100 == 0 && t != 0 )
        {
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }

        // update velocity, v_factor is 'DT / DH / rho'
        // velocity z: vZ = vZ + DT / ( DH * rho ) * A * p;
        vZ += v_factor * A * p;
        // velocity x: vX = vX + DT / ( DH * rho ) * B * p;
        vX += v_factor * B * p;
        // velocity y: vY = vY + DT / ( DH * rho ) * C * p;
        vY += v_factor * C * p;

        // create new Vector(Pointer) with same configuration as vZ
        std::unique_ptr<lama::Vector<ValueType> > helpPtr( vZ.newVector() );
        // get Reference of VectorPointer
        lama::Vector<ValueType>& help = *helpPtr;

        // pressure update
        help =  DH_INV * D * vZ;
        help += DH_INV * E * vX;
        help += DH_INV * F * vY;
        p += p_factor * help; // p_factor is 'DT * M'

        // update seismogram and pressure with source terms
        // CAUTION: elementwise access by setVal and getVal cause performace issues executed on CUDA
        //          should be used rarely
        // TODO: can do this by index operator[] --> no need for DenseVector<>, can use Vector instead
        p.setValue( source_index, p.getValue( source_index ) + source.getValue( t ) );
        seismogram.setValue( t, p.getValue( seismogram_index ) );

        // TODO: plot snapshots of wave propagation ???
    }
}

/*
 *  main for 3-D FD-Algorithm
 *
 *  sets configuration parameters, prints configuration
 *  calculates source
 *  initialize matrices and vectors
 *  do timesteps
 *  writes seismogram data to file
 */
int main( int /*argc*/, char** /*argv[]*/ )
{
    // we do all calculation in double precision

    typedef DefaultReal ValueType;

    // read configuration parameter from file
    Configuration<ValueType> config( "Configuration.txt" );

    // LAMA specific configuration variables

    // execution context
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    hmemo::ContextPtr host = hmemo::Context::getHostPtr();   // explicit host context

    // inter node communicator
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    // inter node distribution
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( config.getN(), comm ) );

    HOST_PRINT( comm, "Acoustic 3-D FD-Algorithm\n\n" );

    if ( comm->getRank() == MASTER )
    {
        config.print();
    }

    // for timing
    double start_t, end_t;

    start_t = common::Walltime::get();
    // get source signal
    // init vector with a sequence of values (MATLAB t=0:DT:(NT*DT-DT);)
    lama::DenseVector<ValueType> source = lama::linearDenseVector<ValueType>( config.getNT(), 0, config.getDT() );
    source.setContextPtr( ctx );
    sourceFunction( source, config.getFC(), config.getAMP(), comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished calculating source in " << end_t - start_t << " sec.\n\n" );

    // printin L2 norm for checking source vector
    //std::cout << "source norm " << source.l2Norm().getValue<ValueType>() << std::endl;

    start_t = common::Walltime::get();
    // calculate sparse matrices
    lama::CSRSparseMatrix<ValueType> A, B, C, D, E, F;
    initializeMatrices( A, B, C, D, E, F, dist, ctx, config.getNX(), config.getNY(), config.getNZ(), comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n" );

    // initialize all needed vectors with zero

    // components of particle velocity

    auto vX = lama::fill<lama::DenseVector<ValueType>>( dist, 0.0, ctx );
    auto vY = lama::fill<lama::DenseVector<ValueType>>( dist, 0.0, ctx );
    auto vZ = lama::fill<lama::DenseVector<ValueType>>( dist, 0.0, ctx );

    // pressure

    auto p = lama::fill<lama::DenseVector<ValueType>>( dist, 0.0, ctx );

    // seismogram data: to store at each time step

    auto seismogram = lama::fill<lama::DenseVector<ValueType>>( config.getNT(), 0.0, host );

    // TODO: load colormap for snapshots ???

    HOST_PRINT( comm, "Start time stepping\n" );

    start_t = common::Walltime::get();
    timesteps( seismogram, source, p, vX, vY, vZ, A, B, C, D, E, F,
               config.getVfactor(), config.getPfactor(), config.getNT(), ValueType( 1 ) / config.getDH(),
               config.getSourceIndex(), config.getSeismogramIndex(), comm, dist );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );

    // print vector data for seismogram plot
    seismogram.writeToFile( "seismogram.mtx" );

    // printing L2 norm for checking seismogram vector
    // std::cout << "L2 Norm Seismogram " << seismogram.l2Norm().getValue<ValueType>() << std::endl;

    return 0;
}

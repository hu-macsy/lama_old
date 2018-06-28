/**
 * @file lama/examples/grid/gridSection.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Example program to work on grid data using reductions
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridSection.hpp>

#include <scai/lama/io/ImageIO.hpp>
#include <scai/lama/io/MatlabIO.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

using namespace common;

typedef double real;
typedef ComplexDouble complex;

SCAI_LOG_DEF_LOGGER( logger, "main" )

#define HOST_PRINT( rank, msg )                 \
    {                                           \
        if ( rank == 0 )                        \
        {                                       \
            std::cout << msg << std::endl;      \
        }                                       \
    }                                           \

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments:
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    // use default communicator, usually MPI

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType np = comm->getSize();    // number of available processors
    IndexType rank = comm->getRank();  // id of this processor, 0 <= rank < np

    if ( argc < 1 )
    {
        HOST_PRINT( rank, "Wrong call, please use : " << argv[0] << " - no arguments" )
    }


    // fixed parameters, do not change otherwise input file Wx.mat does not match

    const IndexType nx = 101; 
    const IndexType nt = 300;
    const real dt = 0.004;
    const real fmin = 3;
    const real fmax = 40;

    // variable parameters, can be modified for benchmarking

    IndexType nz = 240;
    IndexType nsrc = 20;

    //  calculate frequencies
    //  number of positive frequencies and the increment

    const IndexType nf = nt / 2 + 1;
    const real df = 1 / ( nt * dt );

    // define the minimum and maximum frequency number

    IndexType ifmin = static_cast<IndexType>( fmin / df + 0.5 ) + 1;

    if ( ifmin < 2 )
    {
        ifmin = 2;
    }

    IndexType ifmax = static_cast<IndexType>( fmax / df + 0.5 ) + 1;

    if ( ifmax > nf - 1 )
    {
        ifmax = nf - 1;
    }

    IndexType nfeval = ifmax - ifmin + 1;

    HOST_PRINT( rank, "number of frequencies = " << nfeval << ", ifmin = " << ifmin << " - ifmax = " << ifmax )

    GridVector<real> wave( Grid1D( nfeval ), 0 );   // now replicated, will be distributed later

    const real PI = 3.141592653589793;

    // wave is computed on all processors, will be distributed afterwards
    {
        GridWriteAccess<real> wWave( wave );

        IndexType iw = 0;

        for ( IndexType ifr = ifmin; ifr <= ifmax; ++ifr )
        {
            real freq = ( ifr - 1 ) * df;
            real x = ( freq - fmin - df ) / ( 2 * df + fmax - fmin );
            x = common::Math::sin( PI * x );
            wWave[ iw ] = x * x;
            iw++;
        }
    }

    // determine grid sizes, grids with a nfeval dimension are distributed onto available processors

    DistributionPtr distF( new GridDistribution( Grid1D( nfeval ), comm, Grid1D( np ) ) );

    wave.redistribute( distF );  

    DistributionPtr distZXSF( new GridDistribution( Grid4D( nz, nx, nsrc, nfeval ), comm, Grid4D( 1, 1, 1, np ) ) );
    DistributionPtr distXSF( new GridDistribution( Grid3D( nx, nsrc, nfeval ), comm, Grid3D( 1, 1, np ) ) );
    DistributionPtr distXXF( new GridDistribution( Grid3D( nx, nx, nfeval ), comm, Grid3D( 1, 1, np ) ) );

    // So is complex here, avoids later type conversion when added to a complex section

    GridVector<complex> So( distXSF, 0 );   // So( nX, nSrc, nFeval)

    double time = Walltime::get();

    // initialize So:  So( floor( nx / ( nsrc + 1 ) ) * i, i, : ) = wave, i = 1, nsrc

    for ( IndexType isrc = 0; isrc < nsrc; ++isrc )
    {
        IndexType ix = ( nx / ( nsrc + 1 ) ) * ( isrc + 1 ) - 1;

        SCAI_LOG_INFO( logger, "Set So( " << ix << ", " << isrc << ", :) = wave(:)" )

        So( ix, isrc, Range() ) = wave( Range() );
    }

    GridVector<complex> R( Grid2D( nx, nz ), 0 );  // replicated array, each processor has a copy

    SCAI_LOG_INFO( logger, "Reflections at " << ( nz / 2 ) << " and at " << ( nz / 4 ) )

    R( Range(), nz / 2 - 1 ) = 0.2;   // 1st reflection
    R( Range(), nz / 4 - 1 ) = 0.2;   // 2nd reflection

    // load('Wx.mat');  LAMA supports MATLAB files, but only one array, must be set explicitly

    GridVector<complex> Wx( "Wx.mat" );

    // make sure that the read grid data fits our problem here

    SCAI_ASSERT_EQ_ERROR( Wx.globalGrid(), Grid3D( nx, nx, nfeval ), "Input file Wx.mat does not match" )

    // distribute Wx according to the frequencies

    Wx.redistribute( distXXF );
 
    // allocate and initialize arrays Pmin, Pplus

    GridVector<complex> Pmin( distZXSF, 0 );    // shape is ( nz, nx, nsrc, block::nfeval )
    GridVector<complex> Pplus( distZXSF, 0 );   // shape is ( nz, nx, nsrc, block::nfeval )

    // temporary vectors, only declared once

    GridVector<complex> pextr( distXSF, 0 );
    GridVector<complex> pextr_tmp( Grid2D( nx, nsrc ), 0 );
    GridVector<complex> ptmp( Grid2D( nx, nsrc ), 0 );
    GridVector<complex> w_tmp( Grid2D( nx, nx ), 0 );
    GridVector<complex> ri( Grid2D( nx, nx ), 0 );

    // Determine the local range owned by this processor

    IndexType nf_lb = wave.localLB()[0];
    IndexType nfeval_local = wave.localGrid().size( 0 );
    IndexType nf_ub = nf_lb + nfeval_local;

    SCAI_LOG_INFO( logger, *comm << ": my range is : " << nf_lb << ":" << nf_ub << ", nfeval_local = " << nfeval_local << ", nfeval = " << nfeval )
   
    // %% LOOP 1 (can be parallelized in nsrc and nf direction)

    for ( IndexType iz = 0; iz < nz; ++iz )
    {
        SCAI_REGION( "main.loop1" )

        HOST_PRINT( rank, "Loop1 : iter iz = " << iz << " of " << nz )

        SCAI_LOG_INFO( logger, *comm << ": iz = " << iz << " of " << nz )

        for ( IndexType i_f = nf_lb; i_f < nf_ub; i_f++ )
        {
            SCAI_LOG_TRACE( logger, "Loop1: iter ( iz = " << iz << ", i_f = " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            pextr_tmp( Range(), Range() ) = pextr( Range(), Range(), i_f );
            ptmp( Range(), Range() ) = Pmin( iz, Range(), Range(), i_f );

            w_tmp( Range(), Range() ) = Wx( Range(), Range(), i_f );

            // Apply operator r dependent on R (calculations simplified)


            ri.setDiagonal( R( Range(), iz ), 0 );  // ri = diag(R(:,iz));

            ptmp = pextr_tmp - ptmp;
            pextr_tmp.gemm( 1, ri, ptmp );  // pextr_tmp = pextr_tmp + 1 * ri * ptmp

            if ( iz == 0 )
            {
                pextr_tmp( Range(), Range() ) += So( Range(), Range(), i_f );
            }

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );    // ptmp = w_tmp * pextr_tmp;

            pextr( Range(), Range(), i_f ) = ptmp( Range(), Range() );
        }

        Pplus( iz, Range(), Range(), Range() ) = pextr( Range(), Range(), Range() );
    }

    SCAI_LOG_ERROR( logger, *comm << ": Loop 1 -> Loop 2" )

    // %% LOOP 2 (can be parallelized in nsrc and nf direction)

    // pextr=zeros(nx,nsrc,nfeval);

    pextr = 0;
  
    for ( IndexType iz = nz; iz-- > 0;  )
    {
        SCAI_REGION( "main.loop2" )

        HOST_PRINT( rank, "Iter iz = " << iz << " of " << nz )

        for ( IndexType i_f = nf_lb; i_f < nf_ub; ++i_f )
        {
            SCAI_LOG_TRACE( logger, "Loop2: iter ( iz = " << iz << ", i_f = " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            pextr_tmp( Range(), Range() )  = pextr( Range(), Range(), i_f );

            ptmp( Range(), Range() ) = Pplus( iz, Range(), Range(), i_f );

            w_tmp( Range(), Range() ) = Wx( Range(), Range(), i_f );  

            // %Apply operator r dependent on R

            ri.setDiagonal( R( Range(), iz ), 0 );

            // pextr_tmp=pextr_tmp-r*pextr_tmp+r*ptmp;

            ptmp = ptmp - pextr_tmp;
            pextr_tmp.gemm( 1, ri, ptmp );

            // %Apply operator w
            // pextr_tmp=w_tmp*pextr_tmp;

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );

            pextr( Range(), Range(), i_f ) = ptmp( Range(), Range() );

        }  // i_f loop

        Pmin( iz, Range(), Range(), Range() ) = pextr( Range(), Range(), Range() );
    }

    SCAI_LOG_ERROR( logger, *comm << ": Loop 2 -> Loop 3" )

    // %% LOOP 3 (communication necessary after each depth iteration)
    // % Calculate Difference

    // Res=squeeze(Pmin(1,:,:,:));
    // pextr=Res;

    pextr( Range(), Range(), Range() ) = Pmin( 0, Range(), Range(), Range() );

    GridVector<real> grad( Grid2D( nx, nz ), 0 );

    for ( IndexType iz = 0; iz < nz - 1; ++iz )
    {
        SCAI_REGION( "main.loop3" )

        HOST_PRINT( rank, "Iter iz = " << iz << " of " << nz )

        for ( IndexType i_f = nf_lb; i_f < nf_ub; ++i_f )
        {
            SCAI_LOG_TRACE( logger, "Loop3: iter ( iz, i_f ) = ( " << iz << ", " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            pextr_tmp( Range(), Range() ) = pextr( Range(), Range(), i_f );

            // w_tmp=squeeze(Wx(:,:,j));                    %simplified!! w_tmp needs to be calculated for each i 
            // w_tmp = w_tmp'
    
            bool conjFlag = true;
            w_tmp( Range(), Range() ).assignTranspose( Wx( Range(), Range(), i_f ), conjFlag );

            // pextr_tmp=w_tmp'*pextr_tmp;
            // ptmp += w_tmp' * pextr_tmp;

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );

            pextr( Range(), Range(), i_f ) = ptmp( Range(), Range() );
        }

        // %summation over nsrc & nf (problematic???)
        // ptmp=squeeze(Pplus(i+1,:,:,:));
        // grad(:,i+1)=sum(sum(real(pextr.*conj(ptmp)),2),3);

        {
            GridWriteAccess<real> wGrad( grad );
            GridReadAccess<complex> rPextr( pextr );
            GridReadAccess<complex> rPplus( Pplus );

            for ( IndexType ix = 0; ix < nx; ++ix )
            {
                real s = 0;

                for ( IndexType isrc = 0; isrc < nsrc; ++isrc )
                {
                    for ( IndexType i_f = 0; i_f < nfeval_local; ++i_f )
                    {
                        s += common::Math::real( rPextr( ix, isrc, i_f ) * common::Math::conj( rPplus( iz + 1, ix, isrc, i_f ) ) );
                    }
                }

                wGrad( ix, iz + 1 ) = s;
            }
        }
    }
 
    {
        SCAI_REGION( "main.sumReduce" )
        comm->sumArray( grad.getLocalValues() );
    }

    time = Walltime::get() - time;

    std::cout << "Time = " << time << " seconds" << std::endl;

    ImageIO::writeSC( grad, "grad.bmp" );
}

/**
 * @file gridExample.cpp
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

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

using namespace common;

typedef double real;
typedef ComplexDouble complex;

SCAI_LOG_DEF_LOGGER( logger, "main" )

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments:
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( argc < 1 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <inputFileName> <outputFileName>" << std::endl;
    }

    IndexType np = comm->getSize();  // number of available processors

    // fixed parameters

    const IndexType nx = 101;
    const IndexType nt = 300;
    const real dt = 0.004;
    const real fmin = 3;
    const real fmax = 40;

    // variable parameters

    IndexType nz = 120;
    IndexType nsrc = 10;

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

    if ( comm->getRank() == 0 )
    {
    }

    //  calculate sources

    GridVector<real> wave( Grid1D( nfeval ), 0 );   // now replicated, will be distributed later

    const real PI = 3.141592653589793;

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

    DistributionPtr distF( new GridDistribution( Grid1D( nfeval ), comm, Grid1D( np ) ) );

    wave.redistribute( distF );

    // wave.writeToFile( "wave_1.mtx" );

    DistributionPtr distZXSF( new GridDistribution( Grid4D( nz, nx, nsrc, nfeval ), comm, Grid4D( 1, 1, 1, np ) ) );
    DistributionPtr distXSF( new GridDistribution( Grid3D( nx, nsrc, nfeval ), comm, Grid3D( 1, 1, np ) ) );
    DistributionPtr distXXF( new GridDistribution( Grid3D( nx, nx, nfeval ), comm, Grid3D( 1, 1, np ) ) );

    GridVector<real> So( distXSF, 0 );   // So = zeros( nx, nsrc, nfeval);

    // for i = 1 : nsrc
    //    So( floor( nx / ( nsrc + 1 ) ) * i, i, : ) = wave;
    // end

    {   
        // Each processor works on its local parts of So, wave

        GridWriteAccess<real> wSo( So );
        GridReadAccess<real> rWave( wave );

        for ( IndexType isrc = 0; isrc < nsrc; ++isrc )
        {
            IndexType ix = ( nx / ( nsrc + 1 ) ) * ( isrc + 1 ) - 1;

            SCAI_LOG_INFO( logger, "Set So( " << ix << ", " << isrc << ", :) = wave(:)" )

            for ( IndexType i_f = 0; i_f < wSo.size( 2 ); ++i_f )
            {
                wSo( ix, isrc, i_f ) = rWave( i_f );
            }

            // So( ix, isrc, Range() ) = wave( Range() ) = 1.0;
        }
    }

    //  define R
    //  R = zeros( nx, nz );
    //  R( :, floor( nz / 2 ) ) = 0.2;
    //  R( :, floor( nz / 4 ) ) = 0.2;     % add second reflection

    GridVector<real> R( Grid2D( nx, nz ), 0 );  // replicated array, each processor has a copy

    SCAI_LOG_INFO( logger, "Reflection 1 at " << ( nz / 2 ) << ", 2 at " << ( nz / 4 ) )

    {
        GridWriteAccess<real> wR( R );

        IndexType iz1 = nz / 2;
        IndexType iz2 = nz / 4;
 
        for ( IndexType ix = 0; ix < nx; ++ix )
        {
            wR( ix, iz1 - 1 ) = 0.2;
            wR( ix, iz2 - 1 ) = 0.2;
        }
    }

    // be careful: file Wx.mtx contains data in column-major order

    // load Wx matrix
    // MATALB: load('Wx.mat');

    GridVector<complex> Wx( "Wx.mat" );

    // make sure that the read grid data fits our problem here

    SCAI_ASSERT_EQ_ERROR( Wx.globalGrid(), Grid3D( nx, nx, nfeval ), "Input file Wx.mat does not match" )

    // distribute Wx according to the frequencies

    Wx.redistribute( distXXF );
 
    // %% initialise arrays
    // Pmin=zeros(nz,nx,nsrc,nfeval);                     %4D matrices (in general complex)
    // Pplus=zeros(nz,nx,nsrc,nfeval);                    %4D matrices (in general complex)

    GridVector<complex> Pmin( distZXSF, 0 );
    GridVector<complex> Pplus( distZXSF, 0 );

    // temporary vectors, only declared once

    GridVector<complex> pextr( distXSF, 0 );
    GridVector<complex> pextr_tmp( Grid2D( nx, nsrc ), 0 );
    GridVector<complex> ptmp( Grid2D( nx, nsrc ), 0 );
    GridVector<complex> w_tmp( Grid2D( nx, nx ), 0 );
    GridVector<complex> ri( Grid2D( nx, nx ), 0 );

    IndexType nfeval_local = wave.localGrid().size( 0 );

    SCAI_LOG_ERROR( logger, *comm << ": nfeval_local = " << nfeval_local << ", nfeval = " << nfeval )
   
    // %% LOOP 1 (can be parallelized in nsrc and nf direction)

    for ( IndexType iz = 0; iz < nz; ++iz )
    {
        SCAI_LOG_ERROR( logger, *comm << ": iz = " << iz << " of " << nz )

        for ( IndexType i_f = 0; i_f < nfeval_local; i_f++ )
        {
            SCAI_LOG_TRACE( logger, "Loop1: iter ( iz = " << iz << ", i_f = " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            // pextr_tmp=squeeze(pextr(:,:,i_f));

            // ToDo1: GridVector<complex> pextr_tmp( GridSection<complex>( Pextr, Range(), Range(), i_f ) );
            // ToDo2: GridVector<complex> pextr_tmp = Pextr( Range(), Range(), i_f );

            {
                SCAI_LOG_TRACE( logger, "pextr_tmp = pextr" )

                GridReadAccess<complex> rPextr( pextr );
                GridWriteAccess<complex> wPextr_tmp( pextr_tmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr_tmp( ix, isrc ) = rPextr( ix, isrc, i_f  );
                    }
                }
            }

            // ptmp=squeeze(Pmin(iz,:,:,i_f));

            {
                SCAI_LOG_TRACE( logger, "ptmp = Pmin" )

                GridReadAccess<complex> rPmin( Pmin );
                GridWriteAccess<complex> wPtmp( ptmp );

                for ( IndexType ix = 0; ix < nx; ix++ )
                {
                    for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                    {
                        wPtmp( ix, isrc ) = rPmin( iz, ix, isrc, i_f );
                    }
                }
            }

            // w_tmp=squeeze(Wx(:,:,i_f));                  %simplified!! w_tmp needs to be calculated for each i

            {
                SCAI_LOG_TRACE( logger, "w_tmp = Wx(..)" )

                GridReadAccess<complex> rWx( Wx );
                GridWriteAccess<complex> wTmp( w_tmp );

                for ( IndexType ix1 = 0; ix1 < nx; ix1++ )
                {
                    for ( IndexType ix2 = 0; ix2 < nx; ix2++ )
                    {
                        wTmp( ix1, ix2 ) = rWx( ix1, ix2, i_f );
                    }
                }
            }

            // Apply operator r dependent on R (calculations simplified)
            // r=diag(R(:,i));                            %simplified!! in future R can have secondary diagonals or gets close to a full matrix

            {
                SCAI_LOG_TRACE( logger, "ri = diag( R( :, iz = " << iz << " ) )" )

                GridWriteAccess<complex> wRi( ri );
                GridReadAccess<real> rR( R );

                for ( IndexType ix = 0; ix < nx; ++ix )
                {
                    wRi( ix, ix ) = rR( ix, iz );
                }
            }

            ptmp = pextr_tmp - ptmp;
            pextr_tmp.gemm( 1, ri, ptmp );  // pextr_tmp = pextr_tmp + 1 * ri * ptmp

            if ( iz == 0 )
            {
                // sotmp=squeeze(So(:,:,i_f));   % So(nx,nsrc,nfeval)
                // pextr_tmp=pextr_tmp+sotmp;

                {
                    GridReadAccess<real> rSo( So );
                    GridWriteAccess<complex> wPtmp( pextr_tmp );

                    for ( IndexType ix = 0; ix < nx; ++ix )
                    {
                        for ( IndexType isrc = 0; isrc < nsrc; ++isrc )
                        {
                            wPtmp( ix, isrc ) += rSo( ix, isrc, i_f );
                        }
                    }
                }
            }

            // %Apply operator w
            // pextr_tmp=w_tmp*pextr_tmp;

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );


            // pextr(:,:,i_f)=pextr_tmp;
            // pextr( Range(), Range(), i_f ) = ptmp;

            {
                SCAI_LOG_TRACE( logger, "pextr( :, :, " << i_f << " ) = ptmp" )

                GridWriteAccess<complex> wPextr( pextr );
                GridReadAccess<complex> rPtmp( ptmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr( ix, isrc, i_f  ) = rPtmp( ix, isrc );
                    }
                }
            }

        }

        // Pplus( iz, :, :, :) = pextr;  %% pextr( nx, nsrc, nfeval )

        {
            SCAI_LOG_TRACE( logger, "Plus( " << iz << ", :, :, : ) = pextr" )

            GridWriteAccess<complex> wPplus( Pplus );
            GridReadAccess<complex> rPextr( pextr );

            for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
            {
                for ( IndexType ix = 0; ix < nx; ix++ )
                {
                    for ( IndexType i_f = 0; i_f < wPplus.size( 3 ) ; i_f++ )
                    {
                        wPplus( iz, ix, isrc, i_f  ) = rPextr( ix, isrc, i_f );
                    }
                }
            }

            SCAI_LOG_TRACE( logger, "Pplus( " << iz << ", 0, 0, 0 ) = " << wPplus( iz, 0, 0, 0 ) )
        }
    }

    SCAI_LOG_ERROR( logger, *comm << ": Loop 1 -> Loop 2" )

    // %% LOOP 2 (can be parallelized in nsrc and nf direction)

    // pextr=zeros(nx,nsrc,nfeval);

    pextr = 0;
  
    for ( IndexType iz = nz; iz-- > 0;  )
    {
        for ( IndexType i_f = 0; i_f < nfeval_local; ++i_f )
        {
            SCAI_LOG_TRACE( logger, "Loop2: iter ( iz = " << iz << ", i_f = " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            // pextr_tmp = pextr( :, :, i_f );
            {
                GridReadAccess<complex> rPextr( pextr );
                GridWriteAccess<complex> wPextr_tmp( pextr_tmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr_tmp( ix, isrc ) = rPextr( ix, isrc, i_f  );
                    }
                }
            }
           
            // ptmp = squeeze(Pplus(iz,:,:,i_f));

            {
                SCAI_LOG_TRACE( logger, "ptmp = Pplus( " << iz << ", :, :, " << i_f << " )" )

                GridReadAccess<complex> rPplus( Pplus );
                GridWriteAccess<complex> wPtmp( ptmp );

                for ( IndexType ix = 0; ix < nx; ix++ )
                {
                    for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                    {
                        wPtmp( ix, isrc ) = rPplus( iz, ix, isrc, i_f  );
                    }
                }

                SCAI_LOG_TRACE( logger, "Pplus( " << iz << ", 0, 0, 0 ) = " << rPplus( iz, 0, 0, 0 ) )
                SCAI_LOG_TRACE( logger, "ptmp( 0, 0 ) = " << wPtmp( 0, 0 ) )
            }

            // w_tmp=squeeze(Wx(:,:,i_f));                  %simplified!! w_tmp needs to be calculated for each i

            {
                GridReadAccess<complex> rWx( Wx );
                GridWriteAccess<complex> wTmp( w_tmp );

                for ( IndexType ix1 = 0; ix1 < nx; ix1++ )
                {
                    for ( IndexType ix2 = 0; ix2 < nx; ix2++ )
                    {
                        wTmp( ix1, ix2 ) = rWx( ix1, ix2, i_f );
                    }
                }
            }

             // %Apply operator r dependent on R
             // r=diag(R(:,iz));                     %simplified!! in future R can have secondary diagonals or gets close to a full matrix

            {
                GridWriteAccess<complex> wRi( ri );
                GridReadAccess<real> rR( R );

                for ( IndexType ix = 0; ix < nx; ++ix )
                {
                    wRi( ix, ix ) = rR( ix, iz );
                }
            }

            // pextr_tmp=pextr_tmp-r*pextr_tmp+r*ptmp;

            ptmp = ptmp - pextr_tmp;

            pextr_tmp.gemm( 1, ri, ptmp );

            // %Apply operator w
            // pextr_tmp=w_tmp*pextr_tmp;

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );

            // pextr(:,:,i_f)=ptmp;

            {
                GridWriteAccess<complex> wPextr( pextr );
                GridReadAccess<complex> rPtmp( ptmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr( ix, isrc, i_f  ) = rPtmp( ix, isrc );
                    }
                }
            }

        }  // i_f loop

        // Pmin(iz,:,:,:)=pextr;

        {
            GridWriteAccess<complex> wPmin( Pmin );
            GridReadAccess<complex> rPextr( pextr );

            for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
            {
                for ( IndexType ix = 0; ix < nx; ix++ )
                {
                    for ( IndexType i_f = 0; i_f < nfeval_local; i_f++ )
                    {
                        wPmin( iz, ix, isrc, i_f  ) = rPextr( ix, isrc, i_f );
                    }
                }
            }
        }
    }

    SCAI_LOG_ERROR( logger, *comm << ": Loop 2 -> Loop 3" )

    // %% LOOP 3 (communication necessary after each depth iteration)
    // % Calculate Difference

    // Res=squeeze(Pmin(1,:,:,:));
    // pextr=Res;

    {
        GridReadAccess<complex> rPmin( Pmin );
        GridWriteAccess<complex> wPextr( pextr );

        for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
        {
            for ( IndexType ix = 0; ix < nx; ix++ )
            {
                for ( IndexType i_f = 0; i_f < nfeval_local; i_f++ )
                {
                    wPextr( ix, isrc, i_f  ) = rPmin( 0, ix, isrc, i_f );
                }
            }
        }
    }

    GridVector<real> grad( Grid2D( nx, nz ), 0 );

    for ( IndexType iz = 0; iz < nz - 1; ++iz )
    {
        for ( IndexType i_f = 0; i_f < nfeval_local; ++i_f )
        {
            SCAI_LOG_TRACE( logger, "Loop3: iter ( iz, i_f ) = ( " << iz << ", " << i_f 
                                     << " ) of ( " << nz << ", " << nfeval_local << " )" )

            // pextr_tmp=squeeze(pextr(:,:,i_f));

            {
                GridReadAccess<complex> rPextr( pextr );
                GridWriteAccess<complex> wPextr_tmp( pextr_tmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr_tmp( ix, isrc ) = rPextr( ix, isrc, i_f  );
                    }
                }
            }
           
            // w_tmp=squeeze(Wx(:,:,j));                    %simplified!! w_tmp needs to be calculated for each i 
            // w_tmp = w_tmp'

            {
                GridReadAccess<complex> rWx( Wx );
                GridWriteAccess<complex> wTmp( w_tmp );

                for ( IndexType ix1 = 0; ix1 < nx; ix1++ )
                {
                    for ( IndexType ix2 = 0; ix2 < nx; ix2++ )
                    {
                        wTmp( ix2, ix1 ) = common::Math::conj( rWx( ix1, ix2, i_f ) );
                    }
                }
            }

            // pextr_tmp=w_tmp'*pextr_tmp;
            // ptmp += w_tmp' * pextr_tmp;

            ptmp = 0;
            ptmp.gemm( 1, w_tmp, pextr_tmp );

            // pextr(:,:,i_f)=ptmp;

            {
                GridWriteAccess<complex> wPextr( pextr );
                GridReadAccess<complex> rPtmp( ptmp );

                for ( IndexType isrc = 0; isrc < nsrc; isrc++ )
                {
                    for ( IndexType ix = 0; ix < nx; ix++ )
                    {
                        wPextr( ix, isrc, i_f  ) = rPtmp( ix, isrc );
                    }
                }
            }
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

        // sum grad( :, iz + 1 ) over all processors
    }
 
    comm->sumArray( grad.getLocalValues() );

    if ( comm->getRank() == 0 )
    {
        ImageIO::writeSC( grad, "grad.png" );
    }
}
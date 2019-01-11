/**
 * @file examples/bench/fft.cpp
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
 * @brief Benchmark FFT
 * @author Thomas Brandes
 * @date 14.06.2018
 */

#include <scai/lama.hpp>
#include <scai/lama/fft.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

// Matrix & vector related includes

using namespace scai;
using namespace lama;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    SCAI_REGION( "main.main" )

#ifdef SCAI_COMPLEX_SUPPORTED

    typedef DefaultReal ValueType;

    typedef common::Complex<RealType<ValueType>> FFTType;

    // default values for M and N

    IndexType M = 3852;
    IndexType N = 4031;

    if ( argc > 1 )
    {
        M = atoi( argv[1] );
        N = M;
    }

    if ( argc > 2 )
    {
        N = atoi( argv[2] );
    }

    // generate random vector

    hmemo::HArray<ValueType> randomValues( M * N );

    SCAI_REGION_START( "main.setup" )

    double time = common::Walltime::get();

    utilskernel::HArrayUtils::setRandom( randomValues, 1 );

    DenseStorage<ValueType> storage( M, N, randomValues );

    DenseMatrix<ValueType> input( storage );

    // fill up to M2 x N2, dimensions that are the next power 2 values

    const IndexType M2 = 1 << common::Math::nextpow2( M );
    const IndexType N2 = 1 << common::Math::nextpow2( N );

    auto rowDist = std::make_shared<dmemo::BlockDistribution>( M2 );
    auto colDist = std::make_shared<dmemo::BlockDistribution>( N2 );

    // convert to complex matrix, fill it up and distribute it

    auto y = resize<DenseMatrix<FFTType>>( input, rowDist, colDist );

    double setupTime = common::Walltime::get() - time;

    SCAI_REGION_END( "main.setup" )

    std::cout << "Setup: " << setupTime << " seconds" << std::endl;

    time = common::Walltime::get();

    {
        SCAI_REGION( "main.fft" )

        fft( y );        // FFT forward
        ifft( y );       // FFT backward
    }

    double fftTime = common::Walltime::get() - time;

    std::cout << "FFT:   " << fftTime << " seconds" << std::endl;

    SCAI_REGION_START( "main.check" )

    time = common::Walltime::get();

    // resize back to original data

    auto output = resize<DenseMatrix<ValueType>>( y, input.getRowDistributionPtr(), input.getColDistributionPtr() );

    // divide by M * N after fft - ifft to get the original result

    output *= ValueType( 1 ) / ValueType( M2 * N2 );

    auto maxDiff = input.maxDiffNorm( output );

    double checkTime = common::Walltime::get() - time;

    SCAI_REGION_END( "main.check" )

    std::cout << "max diff = " << input.maxDiffNorm( output ) << std::endl;

    std::cout << "Benchmark FFT on 2D matrix " << M << " x " << N << std::endl;
    std::cout << "Setup: " << setupTime << " seconds" << std::endl;
    std::cout << "FFT:   " << fftTime << " seconds" << std::endl;
    std::cout << "Check: " << checkTime << " seconds" << std::endl;

    if ( maxDiff < common::TypeTraits<ValueType>::small() )
    {
        std::cout << "PASSED." << std::endl;
    }
    else
    {
        std::cout << "FAILED." << std::endl;
    }
#else
    std::cout << "FFT benchmark not available, Complex is unsupported." << std::endl;
#endif

}

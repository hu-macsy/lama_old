/**
 * @file fft2D.cpp
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
 * @brief Example program for FFT transformation on pictures
 * @author Thomas Brandes
 * @date 26.03.2018
 */

#include <scai/lama.hpp>

#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/io/ImageIO.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/Settings.hpp>
#include <scai/utilskernel/FFTUtils.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

#ifdef SCAI_COMPLEX_SUPPORTED

   typedef DefaultReal ValueType;
   typedef common::Complex<ValueType> ComplexType;

    if ( argc != 5 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <inputFileName> <outputFileName> lb ub" << std::endl;
        return -1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    IndexType lb = std::stoi( argv[3] );
    IndexType ub = std::stoi( argv[4] );

    // read in the image file, must be a png file

    GridVector<ValueType> image( inputFileName ); 

    std::cout << "read image as grid vector : " << image << std::endl;

    SCAI_ASSERT_EQ_ERROR( image.nDims(), 3, "no color image data" );

    IndexType n1 = image.size( 0 );
    IndexType n2 = image.size( 1 );
    IndexType nc = image.size( 2 );

    SCAI_ASSERT_EQ_ERROR( nc, 3, "only rgb supported" )

    IndexType one = 1;
    IndexType nc1 = one << common::Math::nextpow2( n1 );
    IndexType nc2 = one << common::Math::nextpow2( n2 );

    auto m = zero<DenseStorage<ComplexType>>( nc1, nc2 );

    HArray<ComplexType>& pixels = m.getData();

    std::cout << "Convert image data " << n1 << " x " << n2 
              << " to pixel array " << nc1 << " * " << nc2 << " is " << pixels << std::endl;

    {
        WriteAccess<ComplexType> wPixels( pixels );
        GridReadAccess<ValueType> rImage( image );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                wPixels[ i * nc2 + j ] = rImage( i, j, 0 );
            }
        }
    }

    ValueType maxVal = utilskernel::HArrayUtils::max( image.getLocalValues() );
    ValueType minVal = utilskernel::HArrayUtils::min( image.getLocalValues() );

    std::cout << "pixels range from " << minVal << " - " << maxVal << std::endl << std::endl;

    utilskernel::HArrayUtils::setScalar<ComplexType>( pixels, maxVal, common::BinaryOp::DIVIDE );

    // HArray<ComplexType> tmp;
    // utilskernel::FFTUtils::fft_many( tmp, pixels, nc1, nc2, 1 );
    // utilskernel::HArrayUtils::assign( pixels, tmp );

    int forward = 1;

    utilskernel::FFTUtils::fft_many( pixels, pixels, nc1, nc2, forward );
    m.transposeImpl();
    utilskernel::FFTUtils::fft_many( pixels, pixels, nc2, nc1, forward );
    m.transposeImpl();

    // Operate on FFT data, delete high frequencies

    {
        WriteAccess<ComplexType> wPixels( pixels );

        for ( IndexType i = 0; i < nc1; ++i )
        {
            if ( i >= lb && i <= ub ) continue;

            for ( IndexType j = 0; j < nc2; ++j )
            {
                if ( j >= lb && j <= ub ) continue;

                wPixels[ i * nc2 + j ] = ComplexType( 0 );
            }
        }
    }

    std::cout << "Convert pixel array back to image data " << n1 << " x " << n2 << std::endl;

    int backward = -1;

    utilskernel::FFTUtils::fft_many( pixels, pixels, nc1, nc2, backward );
    m.transposeImpl();
    utilskernel::FFTUtils::fft_many( pixels, pixels, nc2, nc1, backward );
    m.transposeImpl();

    GridVector<ValueType> imageOut( common::Grid3D( n1, n2, 3 ) );

    {
        GridWriteAccess<ValueType> wImage( imageOut );
        ReadAccess<ComplexType> rPixels( pixels );

        for ( IndexType i = 0; i < n1; ++i )
        {
            for ( IndexType j = 0; j < n2; ++j )
            {
                ValueType val = rPixels[ nc2 * i + j ] / static_cast<ValueType>( nc1 * nc2 );
                val *= maxVal;
                wImage( i, j, 0 ) = val;
                wImage( i, j, 1 ) = val;
                wImage( i, j, 2 ) = val;
            }
        }
    }

    imageOut.writeToFile( outputFileName );

#else

    std::cout << "FFT not supported, no Complex type supported" << std::endl;
 
#endif

}

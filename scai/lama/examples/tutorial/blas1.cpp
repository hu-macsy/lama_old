/**
 * @file blas1.cpp
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
 * @brief blas1.cpp
 * @author lschubert
 * @date 17.05.2013
 */

/**
 * NOTE: this example is derived from the blas1 tutorial of ViennaCL
 */

// include necessary system headers
#include <iostream>
#include <time.h>
#include <cstdlib>

// include general lama header
#include <scai/lama.hpp>

//include basic scalar and vector types of LAMA
#include <scai/lama/Scalar.hpp>
#include <scai/lama/DenseVector.hpp>

//include for using different Contexts (e.g. CUDA, default: Host)
#include <scai/hmemo/Context.hpp>

//include for using the NoDistribution
#include <scai/dmemo/NoDistribution.hpp>

//include the generic norm functions of LAMA
#include <scai/lama/norm/all.hpp>

using namespace scai;

int main()
{
    common::Math::srandom( ( unsigned int )time( NULL ) );
    //
    // Define the ValueType used for the vector
    // Change this type definition to double if your gpu supports that
    //
    typedef DefaultReal ScalarType;
    /////////////////////////////////////////////////
    ///////////// Scalar operations /////////////////
    /////////////////////////////////////////////////
    //
    // Define a few scalars:
    //
    ScalarType s1 = 3.1415926;   // static_cast is only needed to switch between t and double by typedef
    ScalarType s2 = 2.71763;
    ScalarType s3 = 42.0;
    // pure scalar operations only can be executed on the host
    std::cout << "Manipulating a few scalars..." << std::endl;
    std::cout << "operator +=" << std::endl;
    s1 += s2;
    std::cout << "operator *=" << std::endl;
    s1 *= s2;
    std::cout << "operator -=" << std::endl;
    s1 -= s2;;
    std::cout << "operator /=" << std::endl;
    s1 /= s2;
    std::cout << "operator +" << std::endl;
    s1 = s2 + s3;
    std::cout << "multiple operators" << std::endl;
    s1 = s2 + s3 * s2 - s3 / s1;
    //
    // Output stream is overloaded as well:
    //
    std::cout << "Scalar s3: " << s3 << std::endl;
    /////////////////////////////////////////////////
    ///////////// Vector operations /////////////////
    /////////////////////////////////////////////////
    //
    // Define a few vectors and fill them with random values
    //
    ScalarType plain_vec[ 10 ];

    for ( unsigned int i = 0; i < 10; ++i )
    {
        plain_vec[ i ] = static_cast<ScalarType>( rand() ) / static_cast<ScalarType>( RAND_MAX );
    }

    lama::DenseVector<ScalarType> lama_vec1;
    lama_vec1.setRawData( 10, plain_vec );
    hmemo::HArray<ScalarType> lama_array1 ( 10, plain_vec );
    auto lama_vec2 = lama::denseVector<ScalarType>( 10, 0 );
    lama_vec2.setDenseValues( lama_array1 );
    lama::DenseVector<ScalarType> lama_vec3( lama_array1 );
    std::cout << "DenseVector with rand values filled" << std::endl;
    //
    // Define the vectors to be used on GPU (CUDA context on device 0) and upload them
    //
    hmemo::ContextPtr cudaContext;

    if ( hmemo::Context::canCreate( common::ContextType::CUDA ) )
    {
        cudaContext = hmemo::Context::getContextPtr( common::ContextType::CUDA, 0 );
    }
    else
    {
        cudaContext = hmemo::Context::getContextPtr( common::ContextType::Host );
    }

    lama_vec1.setContextPtr( cudaContext );
    lama_vec2.setContextPtr( cudaContext );
    lama_vec1.prefetch();
    lama_vec1.prefetch();
    lama_vec3.prefetch( cudaContext );
    std::cout << "vectors copied to CUDA context" << std::endl;
    //
    // Compute the inner product of two GPU vectors and write the result to either CPU or GPU
    //
    s1 = lama_vec1.dotProduct( lama_vec2 );
    std::cout << "dot product calculated" << std::endl;
    //
    // Compute norms:
    //
    s1 = lama_vec1.l1Norm();
    s2 = lama_vec2.l2Norm();
    s3 = lama_vec3.maxNorm();
    std::cout << "norms calculated" << std::endl;
    //
    // Plane rotation of two vectors:
    // Computes (x,y) <- (alpha * x + beta * y, -beta * x + alpha * y)
    //
    ScalarType alpha = 1.1;
    ScalarType beta = 2.3;
    lama_vec1 = alpha * lama_vec1 + beta * lama_vec2;
    lama_vec2 = -beta * lama_vec1 + alpha * lama_vec2;
    std::cout << "plain rotation calculated" << std::endl;
    //
    // Swap the content of two vectors without a temporary vector:
    //
    lama_vec1.swap( lama_vec2 );  //swaps all entries in memory
    std::cout << "vectors swapped" << std::endl;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}

/**
 * @file hmemo/examples/cuda/AliasProblem.cu
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
 * @brief Demo on aliasing problem with LAMA arrays
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/common/macros/assert.hpp>

using namespace scai;
using namespace hmemo;

/* --------------------------------------------------------------------- */

template<typename ValueType>
__global__
void add_kernel( ValueType* array, IndexType n )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;
    ValueType one = 1;

    if ( i < n )
    {
        array[i] += one;
    }
}

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "AliasTest" )

using namespace scai;

typedef hmemo::HArray<double> Array;

void add ( Array& res, const Array& a, const Array& b )
{
    SCAI_ASSERT_LE( a.size(), b.size(), "size mismatch" )
    IndexType n = a.size();
    ContextPtr hostCtx = hmemo::Context::getHostPtr();
    SCAI_LOG_INFO( logger, "res = a + b, n = " << n << ", on " << *hostCtx )
    // Be careful: read accesses should appear before write only access
    hmemo::ReadAccess<double>read1( a, hostCtx );
    hmemo::ReadAccess<double>read2( b, hostCtx );
    hmemo::WriteOnlyAccess<double> write( res, hostCtx, n );
    double* resPtr = write.get();
    const double* aPtr = read1.get();
    const double* bPtr = read2.get();

    for ( IndexType i = 0; i < n; ++i )
    {
        resPtr[i] = aPtr[i] + bPtr[i];
    }
}

void add1 ( Array& a )
{
    ContextPtr gpuCtx = hmemo::Context::getContextPtr( common::ContextType::CUDA );
    int n = a.size();
    SCAI_LOG_INFO( logger, "a = a + 1, n = " << n << ", on " << *gpuCtx )
    hmemo::WriteAccess<double> write( a, gpuCtx );
    double* aPtr = write.get();
    const int blockSize = 256;
    const int nblocks   = ( n + blockSize - 1 ) / blockSize;
    dim3 block( blockSize, 1, 1 );
    dim3 grid( nblocks, 1, 1 );
    add_kernel <<< grid, block>>>( aPtr, n );
}

void printIt( const Array& a )
{
    using namespace std;
    ContextPtr hostCtx = hmemo::Context::getHostPtr();
    int n = a.size();
    hmemo::ReadAccess<double>read( a, hostCtx );
    const double* aPtr = read.get();
    cout << "Array a =";

    for ( IndexType i = 0; i < n; ++i )
    {
        cout << " " << aPtr[i];
    }

    cout << endl;
}

int main()
{
    Array a;
    Array b( 10 , 1.0 );
    Array c( 10 , 2.0 );
    add( a, b, c ); // this is okay
    add( a, a, b ); // this crashed in earlier versions
    // now make sure that we have only a valid copy on GPU
    add1( a );  // done on GPU, invalid values at host
    add( a, a, c );  // might use the old invalid Host values if WriteOnlyAccess is done first
    printIt( a );  // should be 1 + 2 + 1 + 1 + 2 = 7
}


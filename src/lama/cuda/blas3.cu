/**
 * @file blas3.cu
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief blas3.cu
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#include <lama/cuda/CUDAError.hpp>

/**
 * @brief gemm kernel with variable blocksize
 *
 * @param[in] m m
 * @param[in] n n
 * @param[in] k k
 * @param[in] alpha alpha
 * @param[in] A A in row-major order
 * @param[in] lda lower dimension of A
 * @param[in] B B in row-major order
 * @param[in] ldb lower dimension of B
 * @param[in] beta beta
 * @param[out] C C in row-major order
 * @param[in] ldc lower dimension of C
 */
template<typename T>
__global__
void gemm_kernel(
    const int m,
    const int n,
    const int k,
    const T alpha,
    const T* const A,
    const int lda,
    const T* const B,
    const int ldb,
    const T beta,
    T* const C,
    const int ldc )
{
    __shared__ T APart[16 * 16 * 2 * 2];
    T ResPart = 0.0;
    unsigned int Calc = 0;
    unsigned int Load = 256;
    const unsigned int numParts = ( ( k - 1 ) / 16 ) + 1;
    const unsigned int Ci = ( blockIdx.x * 16 + threadIdx.y );
    const unsigned int Cj = ( blockIdx.y * 16 + threadIdx.x );
    //LOAD FIRST PART OF B
    {
        const unsigned char coord = threadIdx.y;
        APart[512 + threadIdx.y * 16 + threadIdx.x] = ( Cj < n && coord < k ) ? B[coord * ldb + Cj] : 0.0;
    }
    //LOAD FIRST PART OF A
    {
        const unsigned char coord = threadIdx.x;
        APart[threadIdx.y * 16 + threadIdx.x] = ( Ci < m && coord < k ) ? A[Ci * lda + coord] : 0.0;
    }
    __syncthreads();

    for( unsigned int ii = 0; ii < numParts; ++ii ) //Loop over Parts of A
    {
        {
            //load next Part of A
            const unsigned char coord = ( ii + 1 ) * 16 + threadIdx.x;
            APart[Load + threadIdx.y * 16 + threadIdx.x] = ( Ci < m && coord < k ) ? A[Ci * lda + coord] : 0.0;
        }
        //#pragma unroll
        //        for(int jj=0; jj < BLOCKSIZEX/BLOCKSIZEY && (jj+ii*BLOCKSIZEX/BLOCKSIZEY)*BLOCKSIZEY<k; ++jj)
        {
            //This loop is not needed for Blocksize: 16x16
            //            if((1+ii)*16<k)
            {
                //load next Part of B if necessary
                const unsigned char coord = ( 1 + ii ) * 16 + threadIdx.y;
                APart[512 + Load + threadIdx.y * 16 + threadIdx.x] =
                    ( Cj < n && coord < k ) ? B[coord * ldb + Cj] : 0.0;//B[Bi*ldb+Bj];
            }
            //            if(Ci<m && Cj<n)
#pragma unroll 16 //Multiplicate Current Tiles
            for( int kk = 15; kk >= 0; --kk )
                //unrolled by hand
            {
                ResPart += APart[threadIdx.y * 16 + Calc + kk] * APart[512 + Calc + kk * 16 + threadIdx.x];
            }

            //            if(ii==numParts-1)
            //                break;
            Load = Calc;
            Calc = 256 - Load;
            __syncthreads();
        }
    }

    {
        //Write Back Result
        if( Ci < m && Cj < n )
        {
            C[Ci * ldc + Cj] = alpha * ResPart + beta * C[Ci * ldc + Cj];
        }
    }
}

template<typename T>
void gemm_launcher(
    const char transA_char,
    const char transB_char,
    const int m,
    const int n,
    const int k,
    const T alpha,
    const T* const A,
    const int lda,
    const T* const B,
    const int ldb,
    const T beta,
    T* const C,
    const int ldc,
    cudaStream_t cs )
{
    LAMA_CHECK_CUDA_ACCESS

    if ( transA_char != transB_char && transA_char != 'N' )
    {
        LAMA_THROWEXCEPTION( "trans = " << transA_char << " not supported for gemm" )
    }

    dim3 dimBlock( 16, 16, 1 );
    dim3 dimGrid( m / ( 16 ) + ( ( m % ( 16 ) ) ? 1 : 0 ), n / 16 + ( ( n % 16 ) ? 1 : 0 ) );
    gemm_kernel <<< dimGrid, dimBlock, 0, cs>>>( m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
}

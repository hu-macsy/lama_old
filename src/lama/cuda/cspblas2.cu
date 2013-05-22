/**
 * @file cspblas2.cu
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
 * @brief cspblas2.cu
 * @author lschubert
 * @date 06.07.2012
 * $Id$
 */

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>

#define IS_ZERO   0x0
#define IS_ONE    0x1
#define IS_MONE   0x2
#define IS_ANY    0x3

/* Define positions for different flags:
 flag 0 : texture
 flag 1 : y equal z
 flag 2, 3 : alpha = 0, 1, -1, else
 flag 4, 5 : beta = 0, 1, -1, else
 */

#define SHIFT_EQUAL  1
#define SHIFT_ALPHA  2
#define SHIFT_BETA   4

#define NO_TEXTURE  0x0
#define USE_TEXTURE 0x1

#define Y_EQ_Z   (0x0 << SHIFT_EQUAL)
#define Y_NE_Z   (0x1 << SHIFT_EQUAL)

#define ALPHA_0  (IS_ZERO << SHIFT_ALPHA)
#define ALPHA_1  (IS_ONE  << SHIFT_ALPHA)
#define ALPHA_M1 (IS_MONE << SHIFT_ALPHA)
#define ALPHA    (IS_ANY  << SHIFT_ALPHA)

#define BETA_0  (IS_ZERO  << SHIFT_BETA)
#define BETA_1  (IS_ONE   << SHIFT_BETA)
#define BETA_M1 (IS_MONE  << SHIFT_BETA)
#define BETA    (IS_ANY   << SHIFT_BETA)

/** Help routine to get flags for a special value. */

template<typename T>
int getSelector( T val )
{
    if ( val == 0 )
    {
        return IS_ZERO;
    }

    if ( val == 1.0 )
    {
        return IS_ONE;
    }

    if ( val == -1.0 )
    {
        return IS_MONE;
    }

    return IS_ANY;
}

/* ========================================================================== */
/*                      CSRSparseMatrix*Vector - KERNEL                       */
/*                                                                            */
/* ========================================================================== */
texture<float,1> texCSRSXref;

texture<int2,1> texCSRDXref;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_CSRx( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_CSRx<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texCSRSXref, i );
}

template<>
__inline__ __device__
double fetch_CSRx<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texCSRDXref, i );
    return __hiloint2double( v.y, v.x );
}

//0. General kernel for
//y_d = alpha * A * x_d + beta * z_d
template<typename T,bool useTexture>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T alpha,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        T summand = beta * z_d[i];
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = alpha * value + summand;
    }
}

//1. Special kernels for
//y_d = alpha * A * x_d + beta * z_d
// with y_d != x_d and y_d != z_d

//1.1 Special kernel for
//y_d = alpha * A * x_d + beta * z_d
// with y_d != x_d and y_d != z_d and alpha not in {0,1,-1} and beta in {0,1,-1}
template<typename T,bool useTexture,int beta>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T alpha,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA<T,beta>( alpha, value, z_d[i] );
    }
}

//1.2 Special kernel for
//y_d = alpha * A * x_d + beta * z_d
// with y_d != x_d and y_d != z_d and alpha in {0,1,-1} and beta not in {0,1,-1}
template<typename T,bool useTexture,int alpha>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA<T,alpha>( beta, z_d[i], value );
    }
}

//1.3 Special kernel for
//y_d = alpha * A * x_d + beta * z_d
// with y_d != x_d and y_d != z_d and alpha in {0,1,-1} and beta in {0,1,-1}
template<typename T,bool useTexture,int alpha,int beta>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA<T,alpha,beta>( value, z_d[i] );
    }
}

//2. Special kernels for
//y_d = alpha * A * x_d + beta * z_d
// with y_d != x_d and y_d == z_d
// => y_d = alpha * A * a_x + beta * y_d

//2.1 Special kernel for
//y_d = alpha * A * x_d + beta * y_d
// with y_d != x_d and y_d == z_d and alpha not in {0,1,-1} and beta not in {0,1,-1}
// => y_d = alpha * A * x_d + beta * y_d
template<typename T,bool useTexture>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T alpha,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA( alpha, value, beta, y_d[i] );
    }
}

//2.2 Special kernel for
//y_d = alpha * A * x_d + beta * y_d
// with y_d != x_d and y_d == z_d and alpha not in {0,1,-1} and beta in {0,1,-1}
template<typename T,bool useTexture,int beta>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T alpha,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA<T,beta>( alpha, value, y_d[i] );
    }
}

//2.3 Special kernel for
//y_d = alpha * A * x_d + beta * y_d
// with y_d != x_d and y_d == z_d and alpha in {0,1,-1} and beta not in {0,1,-1}
template<typename T,bool useTexture,int alpha>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
        }

        y_d[i] = lamaDeviceFMA<T,alpha>( beta, y_d[i], value );
    }
}

//2.4 Special kernel for
//y_d = alpha * A * x_d + beta * y_d
// with y_d != x_d and y_d == z_d and alpha in {0,1,-1} and beta in {0,1,-1}
template<typename T,bool useTexture,int alpha,int beta>
__global__
void cspblas_csrgemv_kernel(
    int n,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        T value = 0.0;

        if ( alpha != 0 )
        {
            const int rowStart = ia_d[i];
            const int rowEnd = ia_d[i + 1];

            for ( int jj = rowStart; jj < rowEnd; ++jj )
            {
                value += a_d[jj] * fetch_CSRx<T,useTexture>( x_d, ja_d[jj] );
            }
        }

        y_d[i] = lamaDeviceFMA<T,alpha,beta>( value, y_d[i] );
    }
}

template<typename T,bool async>
void lama_CSPBLAS_CSRAGEMVPBV_launcher(
    char* transa,
    int* m,
    const T alpha,
    const T* const a_d,
    const int* const ia_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d,
    cudaStream_t stream )
{
    if ( *m == 0 )
    {
        // might happen, handle it here as x_d == y_d would give true as both might be NULL
        return;
    }

    typedef T ValueType;

    if ( x_d == y_d )
    {
        LAMA_THROWEXCEPTION( "INPUTVECTOR_IS_ALIASED_WITH_RESULT_VECTOR" )
    }

    const int block_size = 256;

    dim3 dimBlock( block_size, 1, 1 )

    dim3 dimGrid = makeGrid( *m, dimBlock.x );

    bool useTexture = false; //lama_getUseTex_cuda();

    //TODO:
//    if ( !async && lama_getComputeCapabilityMajor_cuda() < 2 )
//    {
//        //TODO: we need nnc to check if a texture for x_d does not excced the
//        //      max texture size
//        //if(nnc < 134217728)//512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
//        //{
//            useTexture = true;//lama_getUseTex_cuda();
//        //}
//    }
    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texCSRSXref, x_d ), "" );
        }
        else if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texCSRDXref, x_d ), "" );
        }
    }

    LAMA_CUDA_RT_CALL( cudaThreadSetCacheConfig( cudaFuncCachePreferL1 ), "" );
    int selector = 0;

    if ( *transa == 'N' || *transa == 'n' )
    {
        selector = useTexture ? USE_TEXTURE : NO_TEXTURE;
        selector |= ( y_d == z_d ) ? Y_EQ_Z : Y_NE_Z;
        selector |= getSelector<T>( alpha ) << SHIFT_ALPHA;
        selector |= getSelector<T>( beta ) << SHIFT_BETA;
    }
    else if ( *transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c' )
    {
        //TODO: Implement this.
        LAMA_THROWEXCEPTION( "Not implemented yet" )
    }

    // Oops, sorry, but we really have 2 x 2 x 4 x 4 = 64 different cases

    switch ( selector )
    {
    case NO_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 0, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d )
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 0, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, 0, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA:
        cspblas_csrgemv_kernel<ValueType, false, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 1, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 1, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, 1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA:
        cspblas_csrgemv_kernel<ValueType, false, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_0:
        cspblas_csrgemv_kernel < ValueType, false, -1, 0 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_1:
        cspblas_csrgemv_kernel < ValueType, false, -1, 1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, -1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA:
        cspblas_csrgemv_kernel < ValueType, false, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case NO_TEXTURE | Y_EQ_Z | ALPHA | BETA:
        cspblas_csrgemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 0, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 0, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, 0, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_0 | BETA:
        cspblas_csrgemv_kernel<ValueType, false, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 1, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 1, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, 1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_1 | BETA:
        cspblas_csrgemv_kernel<ValueType, false, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_0:
        cspblas_csrgemv_kernel < ValueType, false, -1, 0 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_1:
        cspblas_csrgemv_kernel < ValueType, false, -1, 1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, -1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA:
        cspblas_csrgemv_kernel < ValueType, false, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA | BETA_0:
        cspblas_csrgemv_kernel<ValueType, false, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA | BETA_1:
        cspblas_csrgemv_kernel<ValueType, false, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, false, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case NO_TEXTURE | Y_NE_Z | ALPHA | BETA:
        cspblas_csrgemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
        // now same with USE_TEXTURE
    case USE_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 0, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 0, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, 0, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_0 | BETA:
        cspblas_csrgemv_kernel<ValueType, true, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 1, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 1, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, 1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_1 | BETA:
        cspblas_csrgemv_kernel<ValueType, true, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_0:
        cspblas_csrgemv_kernel < ValueType, true, -1, 0 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_1:
        cspblas_csrgemv_kernel < ValueType, true, -1, 1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, -1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA_M1 | BETA:
        cspblas_csrgemv_kernel < ValueType, true, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, y_d );
        break;
    case USE_TEXTURE | Y_EQ_Z | ALPHA | BETA:
        cspblas_csrgemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, beta, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 0, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 0, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_0 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, 0, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_0 | BETA:
        cspblas_csrgemv_kernel<ValueType, true, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 1, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 1, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, 1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_1 | BETA:
        cspblas_csrgemv_kernel<ValueType, true, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_0:
        cspblas_csrgemv_kernel < ValueType, true, -1, 0 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_1:
        cspblas_csrgemv_kernel < ValueType, true, -1, 1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, -1, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA_M1 | BETA:
        cspblas_csrgemv_kernel < ValueType, true, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA | BETA_0:
        cspblas_csrgemv_kernel<ValueType, true, 0> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA | BETA_1:
        cspblas_csrgemv_kernel<ValueType, true, 1> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA | BETA_M1:
        cspblas_csrgemv_kernel < ValueType, true, -1 > <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, z_d, y_d );
        break;
    case USE_TEXTURE | Y_NE_Z | ALPHA | BETA:
        cspblas_csrgemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>
        ( *m, alpha, a_d, ia_d, ja_d, x_d, beta, z_d, y_d );
        break;
    default:
        THROW_LAMAEXCEPTION( "LAMA_STATUS_NOT_IMPLEMENTED" );
    }

    if ( !async )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "SCSRAGEMVPBV_CUDAKERNEL_FAILED" );

        if ( useTexture )
        {
            if ( sizeof(ValueType) == sizeof(float) )
            {
                LAMA_CUDA_RT_CALL( cudaUnbindTexture( texCSRSXref ), "unbind texture failed" );
            }
            else if ( sizeof(ValueType) == sizeof(double) )
            {
                LAMA_CUDA_RT_CALL( cudaUnbindTexture( texCSRDXref ), "unbind texture failed" );
            }
        }
    }
}

/* ========================================================================== */
/*                      ELLSparseMatrix*Vector - KERNEL                      */
/*                                                                            */
/* ========================================================================== */

texture<float,1> texELLSXref;

texture<int2,1> texELLDXref;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_ELLx( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_ELLx<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texELLSXref, i );
}

template<>
__inline__ __device__
double fetch_ELLx<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texELLDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename T,bool useTexture>
__global__
void cspblas_ellagemvpbv_kernel(
    int n,
    T alpha,
    int nnr,
    const T* a_d,
    const int* ja_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        T summand = 0.0;

        if ( beta != 0.0 )
        {
            summand = beta * z_d[i];
        }

        a_d += i;
        ja_d += i;
        T value = 0.0;

        for ( int kk = 0; kk < nnr; ++kk )
        {
            const T aValue = *a_d;
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += aValue * fetch_ELLx<T,useTexture>( x_d, *ja_d );
            a_d += n;
            ja_d += n;
        }

        y_d[i] = alpha * value + summand;
    }
}

template<typename T,bool async>
void lama_CSPBLAS_ELLAGEMVPBV_launcher(
    char* transa,
    const T alpha,
    int* m,
    int* nnr,
    const T* const a_d,
    const int* const ja_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d,
    cudaStream_t stream )
{
    // TODO read nnc ?? TB
    //    if(nnc>134217728)//512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
    //    {
    //        lama_setLastError(LAMA_STATUS_INPUTVECTOR_EXCEEDS_TEXTURESIZE);
    //        return;
    //    }
    // x_d and y_d need to be not aliased
    // TODO: assert( x_d != y_d );
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( *m, dimBlock.x );
    //TODO: Determine this depending on the compute capability
    bool useTexture = true; //lama_getUseTex_cuda();

    if ( async )
    {
        //useTexture = false;
    }

    if ( useTexture )
    {
        if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, x_d ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }
        else if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLDXref, x_d ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_ellagemvpbv_kernel<T, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_ellagemvpbv_kernel<T, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
    }

    if ( *transa == 'N' || *transa == 'n' )
    {
        if ( useTexture )
        {
            cspblas_ellagemvpbv_kernel<T, true> <<< dimGrid, dimBlock, 0, stream>>>
            ( *m, alpha, *nnr, a_d, ja_d, x_d, beta, z_d, y_d );
        }
        else
        {
            cspblas_ellagemvpbv_kernel<T, false> <<< dimGrid, dimBlock, 0, stream>>>
            ( *m, alpha, *nnr, a_d, ja_d, x_d, beta, z_d, y_d );
        }
    }
    else if ( *transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c' )
    {
//TODO: Implement this.
        throw lama::Exception( lama_getErrorString( LAMA_STATUS_NOT_IMPLEMENTED ) );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" );

    if ( !async )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" );
    }

    {
        if ( useTexture )
        {
            if ( sizeof(T) == sizeof(float) )
            {
                LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
            }
            else if ( sizeof(T) == sizeof(double) )
            {
                LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLDXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
            }
        }
    }
}

template<typename T,bool useTexture>
__global__
void cspblas_ellagemvpbsv_kernel(
    int n,
    T alpha,
    int nnr,
    const int* const ia_d,
    const T* a_d,
    const int* ja_d,
    const int* const rows_d,
    const int nzr,
    const T* const x_d,
    T* y_d )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < nzr )
    {
        const int i = rows_d[id];
        a_d += i;
        ja_d += i;
        T value = 0.0;
        const int noneZeros = ia_d[i];

        for ( int kk = 0; kk < noneZeros; ++kk )
        {
            const T aValue = *a_d;
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += aValue * fetch_ELLx<T,useTexture>( x_d, *ja_d );
            a_d += n;
            ja_d += n;
        }

        y_d[i] += alpha * value;
    }
}

template<typename T>
void lama_CSPBLAS_ELLAGEMVPBSV_launcher(
    char* transa,
    const T alpha,
    int* m,
    int* nnr,
    const int* const ia_d,
    const T* const a_d,
    const int* const ja_d,
    const int* const rows_d,
    const int nzr,
    const T* const x_d,
    T* y_d )
{
//TODO read nnc
//    if(nnc>134217728)//512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
//    {
//        lama_setLastError(LAMA_STATUS_INPUTVECTOR_EXCEEDS_TEXTURESIZE);
//        return;
//    }
//x_d and y_d need to be not aliased
//TODO: assert(x_d != y_d);
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nzr, dimBlock.x );
//TODO: Determine this depending on the compute capability
    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, x_d ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_ellagemvpbv_kernel<T, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_ellagemvpbv_kernel<T, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
    }

    if ( *transa == 'N' || *transa == 'n' )
    {
        if ( useTexture )
        {
            cspblas_ellagemvpbsv_kernel<T, true> <<< dimGrid, dimBlock>>>( *m, alpha, *nnr, ia_d, a_d, ja_d, rows_d, nzr, x_d, y_d );
        }
        else
        {
            cspblas_ellagemvpbsv_kernel<T, false> <<< dimGrid, dimBlock>>>( *m, alpha, *nnr, ia_d, a_d, ja_d, rows_d, nzr, x_d, y_d );
        }
    }
    else if ( *transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c' )
    {
//TODO: Implement this.
        throw lama::Exception( lama_getErrorString( LAMA_STATUS_NOT_IMPLEMENTED ) );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" );

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
    }
}

/* ========================================================================== */
/*                      JDSSparseMatrix*Vector - KERNEL                       */
/*                                                                            */
/* ========================================================================== */
texture<float,1> texJDSSXref;

texture<int2,1> texJDSDXref;

texture<int,1> texJDSdlgRef;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_JDSx( const T* const x, const int i )
{
    return x[i];
}

template<bool useTexture,bool useSharedMemory>
__inline__ __device__
int fetch_JDSdlg( const int* const dlg_d, int[], const int i )
{
    return dlg_d[i];
}

template<>
__inline__ __device__
float fetch_JDSx<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texJDSSXref, i );
}

template<>
__inline__ __device__
double fetch_JDSx<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texJDSDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetch_JDSdlg<true,false>( const int* const, int[], const int i )
{
    return tex1Dfetch( texJDSdlgRef, i );
}

template<>
__inline__ __device__
int fetch_JDSdlg<true,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<>
__inline__ __device__
int fetch_JDSdlg<false,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<typename T,bool useTexture,bool useSharedMem>
__global__
void cspblas_jdsgemv_kernel(
    int n,
    const T alpha,
    const T* const a_d,
    const int* const dlg_d,
    const int ndlg,
    const int* const,
    const int* const ja_d,
    const int* const perm_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* const y_d )
{
    extern __shared__ int dlg[];
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        if ( threadIdx.x < ndlg )
        {
            dlg[threadIdx.x] = dlg_d[threadIdx.x];
        }

        __syncthreads();
    }

    if ( i < n )
    {
        int perm = perm_d[i];
        T summand = 0.0;

        if ( beta != 0.0 )
        {
            summand = beta * z_d[perm];
        }

        T value = 0.0;
        int j;
        int k = i;

        for ( int jj = 0; jj < ndlg; ++jj )
        {
            const int incr = fetch_JDSdlg<useTexture,useSharedMem>( dlg_d, dlg, jj );

            if ( i < incr )
            {
                j = ja_d[k];
                value += a_d[k] * fetch_JDSx<T,useTexture>( x_d, j );
                k += incr;
            }
            else
            {
                break;
            }
        }

//j = ja_d[k];
//value += a_d[k] * fetch_JDSx<T,useTexture>( x_d, j );
        y_d[perm] = alpha * value + summand;
    }
}

template<typename T>
void lama_CSPBLAS_JDSAGEMVPBV_cuda(
    const char* const transa,
    const T alpha,
    const int* const m,
    const T* const a_d,
    const int* const dlg_d,
    const int ndlg,
    const int* const ilg_d,
    const int* const ja_d,
    const int* const perm_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* const y_d )
{
    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( *m, dimBlock.x );
    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSSXref, x_d ), LAMA_STATUS_CUDA_BINDTEX_FAILED );

        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSdlgRef, dlg_d ), LAMA_STATUS_CUDA_BINDTEX_FAILED );
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_jdsgemv_kernel<T, true, false>, cudaFuncCachePreferL1 ),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_jdsgemv_kernel<T, true, true>, cudaFuncCachePreferL1 ),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
    }
    else
    {
        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_jdsgemv_kernel<T, false, false>, cudaFuncCachePreferL1 ),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_jdsgemv_kernel<T, false, true>, cudaFuncCachePreferL1 ),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
    }

    if ( *transa == 'N' || *transa == 'n' )
    {
        if ( useTexture )
        {
            if ( ndlg > block_size )
            {
                cspblas_jdsgemv_kernel<T, true, false> <<< dimGrid, dimBlock>>>( *m, alpha, a_d, dlg_d, ndlg, ilg_d, ja_d,
                        perm_d, x_d, beta, z_d, y_d );
            }
            else
            {
                const int sharedMemSize = ndlg * sizeof(int);
                cspblas_jdsgemv_kernel<T, true, true> <<< dimGrid, dimBlock, sharedMemSize>>>( *m, alpha, a_d, dlg_d, ndlg, ilg_d, ja_d,
                        perm_d, x_d, beta, z_d, y_d );
            }
        }
        else
        {
            if ( ndlg > block_size )
            {
                cspblas_jdsgemv_kernel<T, false, false> <<< dimGrid, dimBlock>>>( *m, alpha, a_d, dlg_d, ndlg, ilg_d, ja_d,
                        perm_d, x_d, beta, z_d, y_d );
            }
            else
            {
                const int sharedMemSize = ndlg * sizeof(int);
                cspblas_jdsgemv_kernel<T, false, true> <<< dimGrid, dimBlock, sharedMemSize>>>( *m, alpha, a_d, dlg_d, ndlg, ilg_d, ja_d,
                        perm_d, x_d, beta, z_d, y_d );
            }
        }
    }
    else if ( *transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c' )
    {
//TODO: Implement this.
        throw lama::Exception( lama_getErrorString( LAMA_STATUS_NOT_IMPLEMENTED ) );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_SJDSAGEMVPBV_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), LAMA_STATUS_SJDSAGEMVPBV_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSSXref ), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );

        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSdlgRef ), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
    }
}

/* ========================================================================== */
/*                      BCSRSparseMatrix*Vector - KERNEL                      */
/*                                                                            */
/* ========================================================================== */

//TODO: Determine Texture Usage based on compute capability
texture<float,1> texBCSRSXref;

texture<int2,1> texBCSRDXref;

template<typename T>
__inline__     __device__ T fetch_BCSRx( const int i )
{
    return tex1Dfetch( texBCSRSXref, i );
}

template<>
__inline__ __device__
double fetch_BCSRx<double>( const int i )
{
    int2 v = tex1Dfetch( texBCSRDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename ValueType,unsigned short blockSizeC>
__global__
void cspblas_bcsrgemv_kernel_rowwise(
    const ValueType alpha,
    const int* const bia_d,
    const int* const bja_d,
    const ValueType* const values_d,
    const int nbu,
    const ValueType beta,
    const ValueType* const z_d,
    ValueType* const result_d )
{
    /* This is the same as                                                    */
    /*      extern __shared__ ValueType product[ ];                           */
    /* But the original version caused the following compile exception:       */
    /*                                                                        */
    /* *** start of expression ***                                            */
    /* [lvalue] variable: name = product, node type: array [0] of double      */
    /* *** end of expression ***                                              */
    /* internal error: check_type_of_variable_node:enk_variable has wrong type*/
    /* 1 catastrophic error detected in the compilation of                    */
    /* "/tmp/tmpxft_00001fc6_00000000-4_lama_CSPBLAS_level2_cuda.cpp1.ii".   */
    /* Compilation aborted.                                                   */
    /* nvcc error   : 'cudafe' died due to signal 6                           */
    /* CMake Error at                                                         */
    /* CMakeFiles/cuda_compile_generated_lama_CSPBLAS_level2_cuda.cu.o.cmake:*/
    /* 232 (message): Error generating file                                   */
    /*                                                                        */
    /* This exception was only thrown, if the kernel was invoked for          */
    /* different value types. If one invocation was commented out, the file   */
    /* was compiled.                                                          */
    /* The solution of this problem is from the nvidea forum:                 */
    /* http://forums.nvidia.com/index.php?showtopic=166905                    */
    /*                                                                        */
    /* It seems like this 'bugg' does not exist in Version 2.3 and 3.1 but    */
    /* only in 3.0. That's at least what they say in the forum.               */
    extern __shared__ char __s_array[];
    ValueType* const product = reinterpret_cast<ValueType*>( __s_array );
    /* ====================================================================== */
    /* STEP 01: figuring out the position within the matrix and the subblock  */
    /*          of the matrix and the CUDA-block.                             */
    /*          Additionally declaring helpful variables.                     */
    /* ====================================================================== */
    /* A) declaration of orientation variables.                               */
    const int localThreadId = localThreadId3D( blockDim, threadIdx );
    const int rowInSubBlock = threadIdx.y;
    const int columnInSubBlock = threadIdx.x;
    /* B) additionally declaration of helpfull variables                      */
    const int blockSizeR = blockDim.y;
    const int numElements = blockSizeR * blockSizeC;
    const int blockRow = ( blockIdx.y * gridDim.x + blockIdx.x ) * blockDim.z + threadIdx.z;
    int blockColumnIndex = bia_d[blockRow] / numElements;
    const int blockColEnd = bia_d[blockRow + 1] / numElements;
    const int posInBlock = rowInSubBlock * blockSizeC + columnInSubBlock;
    /* C) continuing declaration of localizing variables.                     */
    const int row = blockRow * blockSizeR + rowInSubBlock;
    /* ====================================================================== */
    /* STEP 02: set product at the position of the thread to zero.            */
    /* ====================================================================== */
    /* A) set product at the position of the thread to zero, even if the      */
    /*    thread was not used, after all, to have valid results.              */
    product[localThreadId] = 0.0;
    /* B) synchronize threads to make sure, product is completely set to zero.*/
    __syncthreads();

    /* C) thread leaves the block, if its block-row is outside of the matrix. */
    if ( blockRow >= nbu )
    {
        return;
    }

    /* ====================================================================== */
    /* STEP 03: multiply in a loop over the blocks of this block column every */
    /*          value with the belonging value of the vector and add it in    */
    /*          the same step to the value of the shared memory.              */
    /* ====================================================================== */
    for ( int i = 0; blockColumnIndex < blockColEnd; ++blockColumnIndex, ++i )
    {
        /* A) declaration of orientation variables.                           */
//        int column      = bja_d[blockColumnIndex]*blockSizeC + columnInSubBlock;
        int column = bja_d[blockColumnIndex] + columnInSubBlock;
        int position = bia_d[blockRow] + i * numElements + posInBlock;
        /* B) multiplication                                                  */
        product[localThreadId] += values_d[position] * fetch_BCSRx<ValueType>( column );
    }

    /* C) synchronizing threads.                                              */
    __syncthreads();

    /* ====================================================================== */
    /* STEP 04: add in a loop over the subblock the values of the subblock.   */
    /* ====================================================================== */

    switch ( blockSizeC )
    {
    case 16:

        if ( columnInSubBlock < 8 )
        {
            product[localThreadId] += product[localThreadId + 8];
        }

        __syncthreads();
        /* intended fall through                                          */
    case 8:

        if ( columnInSubBlock < 4 )
        {
            product[localThreadId] += product[localThreadId + 4];
        }

        __syncthreads();
        /* intended fall through                                          */
    case 4:

        if ( columnInSubBlock < 2 )
        {
            product[localThreadId] += product[localThreadId + 2];
        }

        __syncthreads();
        /* intended fall through                                          */
    case 2:

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 1];
        }

        __syncthreads();
        /* intended fall through                                          */
    case 1:
        /* nothing to do, but avoid stucking within the default-case,     */
        /* just a break.                                                  */
        break;
    case 7:

        if ( columnInSubBlock < 3 )
        {
            product[localThreadId] += product[localThreadId + 4];
        }

        __syncthreads();

        if ( columnInSubBlock < 2 )
        {
            product[localThreadId] += product[localThreadId + 2];
        }

        __syncthreads();

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 1];
        }

        __syncthreads();
        break;
    case 6:

        if ( columnInSubBlock < 3 )
        {
            product[localThreadId] += product[localThreadId + 3];
        }

        __syncthreads();
        /* intended fall through                                          */
    case 3:

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 2];
        }

        __syncthreads();

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 1];
        }

        __syncthreads();
        break;
    case 5:

        if ( columnInSubBlock < 2 )
        {
            product[localThreadId] += product[localThreadId + 3];
        }

        __syncthreads();

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 2];
        }

        __syncthreads();

        if ( columnInSubBlock < 1 )
        {
            product[localThreadId] += product[localThreadId + 1];
        }

        __syncthreads();
        break;
    default:

        for ( int i = blockDim.y; i > 1; i = halve( i ) )
        {
            /* A) check if the subblock-column of the thread is within the*/
            /*    first halve of the interval.                            */
            if ( columnInSubBlock < div2( i ) )
            {
                /* B) add the product-sum of the other thread with the    */
                /*    own.                                                */
                product[localThreadId] += product[localThreadId + halve( i )];
            }

            /* C) synchronize threads.                                    */
            __syncthreads();
        }
    }

    /* ====================================================================== */
    /* STEP 05: threads having column 0 write the sum to global memory.       */
    /* ====================================================================== */

    if ( columnInSubBlock == 0 )
    {
        result_d[row] = alpha * product[localThreadId] + beta * z_d[row];
    }
}

template<typename T>
void lama_CSPBLAS_BCSRAGEMVPBV_cuda(
    const char*,
    const T alpha,
    const int* const bia_d,
    const int* const bja_d,
    const T* const values_d,
    const int blockR,
    const int blockC,
    const int nbu,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* const y_d )
{
    const int numThreadsPerBlock = 128;
    const int numElements = blockR * blockC;
    /* min and max are functions included from cuda_runtime.h                 */
    const int numBlockRowsPerCUDABlock = max( min( numThreadsPerBlock / numElements, 64 ), 1 );
    dim3 dimBlock( blockC, blockR, numBlockRowsPerCUDABlock );
    dim3 dimGrid = makeGrid( nbu * numElements, dimBlock.x * dimBlock.y * dimBlock.z );
    cudaError error = cudaBindTexture( NULL, texBCSRSXref, x_d );

    if ( cudaSuccess != error )
    {
        LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_BINDTEX_FAILED );
        return;
    }

    switch ( blockC )
    {
    case 1:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,1>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }

        cspblas_bcsrgemv_kernel_rowwise<T, 1>
        <<< dimGrid, dimBlock, dimBlock.x * dimBlock.y * dimBlock.z * sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 2:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,2>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 2>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 3:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,3>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 3>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 4:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,4>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 4>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 5:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,5>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 5>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 6:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,6>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 6>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 7:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,7>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 7>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 8:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,8>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 8>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    case 16:
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,16>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 16>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    default:
        /* Of course, 255 is not the real block size. Actually this does  */
        /* not matter, because the kernel is optimized for the block sizes*/
        /* listed above and will run a loop for any other size. But still */
        /* a default-value is needed at compiletime, so 255 was chosen,   */
        /* because it is on most platforms the biggest short value and    */
        /* represents -1.                                                 */
        error = cudaFuncSetCacheConfig( cspblas_bcsrgemv_kernel_rowwise<T,255>, cudaFuncCachePreferShared );

        if ( cudaSuccess != error )
        {
            LAMA_CUDA_RT_CALL( error, LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
            return;
        }
        cspblas_bcsrgemv_kernel_rowwise<T, 255>
        <<< dimGrid, dimBlock, dimBlock.x* dimBlock.y* dimBlock.z* sizeof( T )>>>(
            alpha, bia_d, bja_d, values_d, nbu, beta, z_d, y_d );
        break;
    }

    error = cudaGetLastError();

    if ( cudaSuccess != error )
    {
        LAMA_CUDA_RT_CALL( error, LAMA_STATUS_SBCSRGEMV_CUDAKERNEL_FAILED );
        return;
    }

    error = cudaStreamSynchronize( 0 );

    if ( cudaSuccess != error )
    {
        LAMA_CUDA_RT_CALL( error, LAMA_STATUS_SBCSRGEMV_CUDAKERNEL_FAILED );
        return;
    }

    throw lama::Exception( lama_getErrorString( LAMA_STATUS_SUCCESS ) );
}

/* ========================================================================== */
/*                      DIASparseMatrix*Vector - KERNEL                       */
/*                                                                            */
/* ========================================================================== */
//TODO: Determine Texture Usage based on compute capability
texture<float,1> texDIASXref;

texture<int2,1> texDIADXref;

//DIASparseMatrix
template<typename T>
__inline__     __device__ T fetch_DIAx( const int i )
{
    return tex1Dfetch( texDIASXref, i );
}

template<>
__inline__ __device__
double fetch_DIAx<double>( const int i )
{
    int2 v = tex1Dfetch( texDIADXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename ValueType>
__global__ void cspblas_diaagemvpbv_kernel(
    const int nnu,
    const ValueType alpha,
    const int nnc,
    const int nd,
    const int nnd,
    const int* const ia_d,
    const ValueType* const data_d,
    ValueType beta,
    const ValueType* const z_d,
    ValueType* const result_d,
    const bool first )
{
    extern __shared__ int ia[];

    if ( threadIdx.x < nd )
    {
        ia[threadIdx.x] = ia_d[threadIdx.x];
    }

    __syncthreads();
    int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < nnu )
    {
        if ( first )
        {
            if ( beta != 0.0 )
            {
                result_d[i] = beta * z_d[i];
            }
            else
            {
                result_d[i] = 0.0;
            }
        }

        ValueType temp = 0.0;

        for ( int jj = 0; jj < nd; jj++ )
        {
            int j = i + ia[jj];

            if ( j >= 0 && j < nnc )
            {
                ValueType val = alpha * data_d[nnd * jj + i];
                temp += val * fetch_DIAx<ValueType>( j );
            }
        }

        result_d[i] += temp;
    }
}

template<typename ValueType>
__global__ void cspblas_diaagemvpbv_kernel_small(
    const int nnu,
    const ValueType alpha,
    const int nnc,
    const int nd,
    const int nnd,
    const int* const ia_d,
    const ValueType* const data_d,
    ValueType beta,
    const ValueType* const z_d,
    ValueType* const result_d )
{
    extern __shared__ int ia[];

    if ( threadIdx.x < nd )
    {
        ia[threadIdx.x] = ia_d[threadIdx.x];
    }

    __syncthreads();
    int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < nnu )
    {
        ValueType summand = 0.0;

        if ( beta != 0.0 )
        {
            summand = beta * z_d[i];
        }

        ValueType temp = 0.0;

        for ( int jj = 0; jj < nd; jj++ )
        {
            int j = i + ia[jj];

            if ( j >= 0 && j < nnc )
            {
                ValueType val = data_d[nnd * jj + i];
                temp += val * fetch_DIAx<ValueType>( j );
            }
        }

        result_d[i] = alpha * temp + summand;
    }
}

template<typename T>
void lama_CSPBLAS_DIAAGEMVPBV_cuda(
    const char* transa,
    const T alpha,
    const int* nnu,
    const int* nnc,
    const int* nd,
    const int* nnd,
    const int* const ia_d,
    const T* const data_d,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d )
{
    if ( *nnc > 134217728 ) //512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
    {
        throw lama::Exception( lama_getErrorString( LAMA_STATUS_INPUTVECTOR_EXCEEDS_TEXTURESIZE ) );
    }

    const int block_size = 256;

    dim3 dimBlock( block_size, 1, 1 );

    dim3 dimGrid = makeGrid( *nnu, dimBlock.x );

    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texDIASXref, x_d ), LAMA_STATUS_CUDA_BINDTEX_FAILED );

    LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_diaagemvpbv_kernel<T>, cudaFuncCachePreferShared ),
                       LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );

    LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( cspblas_diaagemvpbv_kernel_small<T>, cudaFuncCachePreferShared ),
                       LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );

    if ( *transa == 'n' || *transa == 'N' )
    {
        if ( block_size <= *nd )
        {
            for ( unsigned int base = 0; base < *nd; base += block_size )
            {
                int num_diags = min( *nd - base, (unsigned int) block_size );
                cspblas_diaagemvpbv_kernel <<< dimGrid, dimBlock, num_diags * sizeof( int )>>>( *nnu, alpha, *nnc,
                        num_diags, *nnd, ia_d + base, data_d + ( base*( *nnu ) ), beta, z_d, y_d, ( base == 0 ) );
            }
        }
        else
        {
            cspblas_diaagemvpbv_kernel_small <<< dimGrid, dimBlock, block_size * sizeof( int )>>>( *nnu, alpha, *nnc,
                    *nd, *nnd, ia_d, data_d, beta, z_d, y_d );
        }
    }
    else if ( *transa == 't' || *transa == 'T' || *transa == 'c' || *transa == 'C' )
    {
        throw lama::Exception( lama_getErrorString( LAMA_STATUS_NOT_IMPLEMENTED ) );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_SDIAAGEMVPBV_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), LAMA_STATUS_SDIAAGEMVPBV_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texDIASXref ), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
}

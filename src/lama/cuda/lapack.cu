/**
 * @file lapack.cu
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
 * @brief lapack.cu
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>

texture<float,1> texCSRJacobiSXref;

texture<int2,1> texCSRJacobiDXref;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_CSRJacobix( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_CSRJacobix<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texCSRJacobiSXref, i );
}

template<>
__inline__ __device__
double fetch_CSRJacobix<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texCSRJacobiDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename T,bool useTexture>
__global__
void cspblas_CSRJacobi_kernel(
    const int* const ia_d,
    const int* const ja_d,
    const T* const a_d,
    const int nnu,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( i < nnu )
    {
        T temp = rhs_d[i];
        const int rowStart = ia_d[i];
        const int rowEnd = ia_d[i + 1];
        const T a_diag = a_d[rowStart];
        for ( int jj = rowStart + 1; jj < rowEnd; ++jj )
        {
            temp -= a_d[jj] * fetch_CSRJacobix<T,useTexture>( u_old_d, ja_d[jj] );
        }
        if ( omega == 0.5 )
        {
            u_d[i] = omega * ( fetch_CSRJacobix<T,useTexture>( u_old_d, i ) + temp / a_diag );
        }
        else if ( omega == 1.0 )
        {
            u_d[i] = temp / a_diag;
        }
        else
        {
            u_d[i] = omega * ( temp / a_diag ) + ( 1.0 - omega ) * fetch_CSRJacobix<T,useTexture>( u_old_d, i );
        }
    }
}

template<typename T>
void CSRJacobiAsync_launcher(
    const int* const ia_d,
    const int* const ja_d,
    const T* const a_d,
    const int nnu,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega,
    cudaStream_t stream )
{
    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nnu, dimBlock.x );

    //TODO: Determine this depending on the comput capability
    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texCSRJacobiSXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_CSRJacobi_kernel<T,true>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_CSRJacobi_kernel<T,false>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }

    if ( useTexture )
    {
        cspblas_CSRJacobi_kernel<T,true><<<dimGrid,dimBlock,0,stream>>>(ia_d,ja_d,a_d,nnu,rhs_d,u_d,u_old_d,omega);
    }
    else
    {
        cspblas_CSRJacobi_kernel<T,false><<<dimGrid,dimBlock,0,stream>>>(ia_d,ja_d,a_d,nnu,rhs_d,u_d,u_old_d,omega);
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_SCSRJACOBI_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiSXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
    }
}

template<typename T,bool useTexture>
__global__
void cspblas_CSRJacobi_halo_kernel(
    T* const u_d,
    const int* const ia_local_d,
    const T* const a_local_d,
    const int* const ia_halo_d,
    const int* const ja_halo_d,
    const T* const a_halo_d,
    const int* const rows_d,
    const int nzr,
    const int n,
    const T* const u_old_d,
    const T omega )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( id < nzr )
    {
        int i = id;
        if ( nzr < n )
        {
            i = rows_d[id];
        }

        T temp = 0.0;
        const int rowStart = ia_halo_d[i];
        const int rowEnd = ia_halo_d[i + 1];

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            temp += a_halo_d[jj] * fetch_CSRJacobix<T,useTexture>( u_old_d, ja_halo_d[jj] );
        }

        const T a_diag = a_local_d[ia_local_d[i]];
        u_d[i] -= temp * omega / a_diag;
    }
}

template<typename T>
void CSRJacobiHalo_launcher(
    T* const u_d,
    const int* const ia_local_d,
    const T* const a_local_d,
    const int* const ia_halo_d,
    const int* const ja_halo_d,
    const T* const a_halo_d,
    const int* const rows_d,
    const int nzr,
    const int n,
    const T* const u_old_d,
    const T omega )
{
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nzr, dimBlock.x );

    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {

        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texCSRJacobiDXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texCSRJacobiSXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_CSRJacobi_halo_kernel<T,true>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_CSRJacobi_halo_kernel<T,false>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }

    if ( useTexture )
    {
        cspblas_CSRJacobi_halo_kernel<T,true><<<dimGrid,dimBlock>>>(u_d,ia_local_d,a_local_d,ia_halo_d,ja_halo_d,a_halo_d,rows_d,nzr,n,u_old_d,omega);
    }
    else
    {
        cspblas_CSRJacobi_halo_kernel<T,false><<<dimGrid,dimBlock>>>(u_d,ia_local_d,a_local_d,ia_halo_d,ja_halo_d,a_halo_d,rows_d,nzr,n,u_old_d,omega);
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
    }
}

texture<float,1> texELLJacobiSXref;

texture<int2,1> texELLJacobiDXref;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_ELLJacobix( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_ELLJacobix<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texELLJacobiSXref, i );
}

template<>
__inline__ __device__
double fetch_ELLJacobix<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texELLJacobiDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename T,bool useTexture>
__global__
void cspblas_ELLJacobi_kernel(
    const int nnr,
    const int* ja_d,
    const T* a_d,
    const int n,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( i < n )
    {
        T temp = rhs_d[i];
        a_d += i;
        ja_d += i;

        const T a_diag = *a_d;
        a_d += n;
        ja_d += n;
        for ( int kk = 1; kk < nnr; ++kk )
        {
            const T aValue = *a_d;
            temp -= aValue * fetch_ELLJacobix<T,useTexture>( u_old_d, *ja_d );
            a_d += n;
            ja_d += n;
        }
        if ( omega == 0.5 )
        {
            u_d[i] = omega * ( fetch_ELLJacobix<T,useTexture>( u_old_d, i ) + temp / a_diag );
        }
        else if ( omega == 1.0 )
        {
            u_d[i] = temp / a_diag;
        }
        else
        {
            u_d[i] = omega * ( temp / a_diag ) + ( 1.0 - omega ) * fetch_ELLJacobix<T,useTexture>( u_old_d, i );
        }
    }
}

template<typename T>
void ELLJacobiAsync_launcher(
    const int nnr,
    const int* const ja_d,
    const T* const a_d,
    const int nnu,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega,
    cudaStream_t stream )
{
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nnu, dimBlock.x );

//TODO: Determine this depending on the comput capability
    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texELLJacobiSXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_ELLJacobi_kernel<T,true>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_ELLJacobi_kernel<T,false>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }

    if ( useTexture )
    {
        cspblas_ELLJacobi_kernel<T,true><<<dimGrid,dimBlock,0,stream>>>(nnr,ja_d,a_d,nnu,rhs_d,u_d,u_old_d,omega);
    }
    else
    {
        cspblas_ELLJacobi_kernel<T,false><<<dimGrid,dimBlock,0,stream>>>(nnr,ja_d,a_d,nnu,rhs_d,u_d,u_old_d,omega);
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_DCSRJACOBI_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture(texELLJacobiSXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
    }
}

template<typename T,bool useTexture>
__global__
void cspblas_ELLJacobi_halo_kernel(
    T* const u_d,
    const T* const a_local_d,
    const int* const ia_halo_d,
    const int* const ja_halo_d,
    const T* const a_halo_d,
    const int* const rows_d,
    const int nzr,
    const int n,
    const T* const u_old_d,
    const T omega )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( id < nzr )
    {
        int i = id;
        if ( nzr < n )
        {
            i = rows_d[id];
        }

        T temp = 0.0;

        int pos = i;
        const int rowEnd = ia_halo_d[i];
        for ( int jj = 0; jj < rowEnd; ++jj )
        {
            temp += a_halo_d[pos] * fetch_ELLJacobix<T,useTexture>( u_old_d, ja_halo_d[pos] );
            pos += n;
        }

        const T a_diag = a_local_d[i];
        u_d[i] -= temp * omega / a_diag;
    }
}

template<typename T>
void ELLJacobiHalo_launcher(
    T* const u_d,
    const T* const a_local_d,
    const int* const ia_halo_d,
    const int* const ja_halo_d,
    const T* const a_halo_d,
    const int* const rows_d,
    const int nzr,
    const int n,
    const T* const u_old_d,
    const T omega )
{
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nzr, dimBlock.x );

    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {

        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texELLJacobiDXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texELLJacobiSXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_ELLJacobi_halo_kernel<T,true>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_ELLJacobi_halo_kernel<T,false>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );

    }

    if ( useTexture )
    {
        cspblas_ELLJacobi_halo_kernel<T,true><<<dimGrid,dimBlock>>>(u_d,a_local_d,ia_halo_d,ja_halo_d,a_halo_d,rows_d,nzr,n,u_old_d,omega);
    }
    else
    {
        cspblas_ELLJacobi_halo_kernel<T,false><<<dimGrid,dimBlock>>>(u_d,a_local_d,ia_halo_d,ja_halo_d,a_halo_d,rows_d,nzr,n,u_old_d,omega);
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
    }
}

texture<float,1> texJDSJacobiSXref;

texture<int2,1> texJDSJacobiDXref;

texture<int,1> texJDSJacobidlgRef;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_JDSJacobix( const T* const x, const int i )
{
    return x[i];
}

template<bool useTexture,bool useSharedMemory>
__inline__ __device__
int fetch_JDSJacobidlg( const int* const dlg_d, int[], const int i )
{
    return dlg_d[i];
}

template<>
__inline__ __device__
float fetch_JDSJacobix<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texJDSJacobiSXref, i );
}

template<>
__inline__ __device__
double fetch_JDSJacobix<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texJDSJacobiDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetch_JDSJacobidlg<true,false>( const int* const, int[], const int i )
{
    return tex1Dfetch( texJDSJacobidlgRef, i );
}

template<>
__inline__ __device__
int fetch_JDSJacobidlg<true,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<>
__inline__ __device__
int fetch_JDSJacobidlg<false,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<typename T,bool useTexture,bool useSharedMem>
__global__
void cspblas_jdsjacobi_kernel(
    const T* const a_d,
    const int* const dlg_d,
    const int ndlg,
    const int* const ja_d,
    const int* const perm_d,
    const int nnu,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega )
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

    if ( i < nnu )
    {
        const int perm = perm_d[i];

        T temp = rhs_d[perm];

        const T a_diag = a_d[i];

        int k = i + fetch_JDSJacobidlg<useTexture,useSharedMem>( dlg_d, dlg, 0 );
        ;
        for ( int jj = 1; jj < ndlg; ++jj )
        {
            const int incr = fetch_JDSJacobidlg<useTexture,useSharedMem>( dlg_d, dlg, jj );
            if ( i < incr )
            {
                const int j = ja_d[k];
                temp -= a_d[k] * fetch_JDSJacobix<T,useTexture>( u_old_d, j );
                k += incr;
            }
            else
            {
                break;
            }
        }

        if ( omega == 0.5 )
        {
            u_d[perm] = omega * ( fetch_JDSJacobix<T,useTexture>( u_old_d, perm ) + temp / a_diag );
        }
        else if ( omega == 1.0 )
        {
            u_d[perm] = temp / a_diag;
        }
        else
        {
            u_d[perm] = omega * ( temp / a_diag ) + ( 1.0 - omega ) * fetch_JDSJacobix<T,useTexture>( u_old_d, perm );
        }
    }
}

template<typename T>
void JDSJacobiAsync_launcher(
    const T* const a_d,
    const int* const dlg_d,
    const int ndlg,
    const int* const /*ilg_d*/,
    const int* const ja_d,
    const int* const perm_d,
    const int nnu,
    const T* const rhs_d,
    T* const u_d,
    const T* const u_old_d,
    const T omega,
    cudaStream_t stream )
{
    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nnu, dimBlock.x );

    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texJDSJacobiSXref,u_old_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );

        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texJDSJacobidlgRef,dlg_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_jdsjacobi_kernel<T,true,false>,cudaFuncCachePreferL1),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_jdsjacobi_kernel<T,true,true>,cudaFuncCachePreferL1),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }

    }
    else
    {
        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_jdsjacobi_kernel<T,false,false>,cudaFuncCachePreferL1),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_jdsjacobi_kernel<T,false,true>,cudaFuncCachePreferL1),
                               LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
        }
    }

    if ( useTexture )
    {
        if ( ndlg > block_size )
        {
            cspblas_jdsjacobi_kernel<T,true,false><<<dimGrid,dimBlock,0,stream>>>(a_d,dlg_d,ndlg,ja_d,perm_d,nnu,rhs_d,u_d,u_old_d,omega);
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof(int);
            cspblas_jdsjacobi_kernel<T,true,true><<<dimGrid,dimBlock,sharedMemSize,stream>>>(a_d,dlg_d,ndlg,ja_d,perm_d,nnu,rhs_d,u_d,u_old_d,omega);
        }
    }
    else
    {
        if ( ndlg > block_size )
        {
            cspblas_jdsjacobi_kernel<T,false,false><<<dimGrid,dimBlock,0,stream>>>(a_d,dlg_d,ndlg,ja_d,perm_d,nnu,rhs_d,u_d,u_old_d,omega);
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof(int);
            cspblas_jdsjacobi_kernel<T,false,true><<<dimGrid,dimBlock,sharedMemSize,stream>>>(a_d,dlg_d,ndlg,ja_d,perm_d,nnu,rhs_d,u_d,u_old_d,omega);
        }
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_SJDSJACOBI_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture(texJDSJacobiSXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        if ( ndlg > block_size )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texJDSJacobidlgRef), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
    }
}

template<typename T,bool useTexture>
__global__
void cspblas_JDSJacobi_halo_kernel(
    const T* const a_local_d,
    const T* const a_halo_d,
    const int* const dlg_halo_d,
    const int ndlg_halo,
    const int* const ilg_halo_d,
    const int* const ja_halo_d,
    const int* const perm_halo_d,
    const int* const rows_d,
    const int nnu,
    T* const u_local_d,
    const T* const u_old_halo_d,
    const T omega )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( id < nnu )
    {
        int i = id;
        if ( rows_d != NULL )
        {
            i = rows_d[id];
        }
        T temp = 0.0;

        int pos = i;
        const int rowEnd = ilg_halo_d[i];
        for ( int jj = 0; jj < rowEnd; ++jj )
        {
            temp += a_halo_d[pos] * fetch_JDSJacobix<T,useTexture>( u_old_halo_d, ja_halo_d[pos] );
            pos += dlg_halo_d[jj];
        }

        const T a_diag = a_local_d[i];
        u_local_d[perm_halo_d[i]] -= temp * omega / a_diag;
    }
}

template<typename T>
void JDSJacobiHalo_launcher(
    const T* const a_local_d,
    const T* const a_halo_d,
    const int* const dlg_halo_d,
    const int ndlg_halo,
    const int* const ilg_halo_d,
    const int* const ja_halo_d,
    const int* const perm_halo_d,
    const int* const rows_d,
    const int nnu,
    T* const u_local_d,
    const T* const u_old_halo_d,
    const T omega )
{
    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( nnu, dimBlock.x );

    const bool useTexture = true; //lama_getUseTex_cuda();

    if ( useTexture )
    {

        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texJDSJacobiDXref,u_old_halo_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture(NULL,texJDSJacobiSXref,u_old_halo_d), LAMA_STATUS_CUDA_BINDTEX_FAILED );
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_JDSJacobi_halo_kernel<T,true>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig(cspblas_JDSJacobi_halo_kernel<T,false>, cudaFuncCachePreferL1),
                           LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED );

    }

    if ( useTexture )
    {
        cspblas_JDSJacobi_halo_kernel<T,true><<<dimGrid,dimBlock>>>(a_local_d,a_halo_d,dlg_halo_d,
                ndlg_halo,ilg_halo_d,ja_halo_d,
                perm_halo_d,rows_d,nnu,
                u_local_d,u_old_halo_d,omega);
    }
    else
    {
        cspblas_JDSJacobi_halo_kernel<T,false><<<dimGrid,dimBlock>>>(a_local_d,a_halo_d,dlg_halo_d,
                ndlg_halo,ilg_halo_d,ja_halo_d,
                perm_halo_d,rows_d,nnu,
                u_local_d,u_old_halo_d,omega);
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED );

    if ( useTexture )
    {
        if ( sizeof(T) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
        else if ( sizeof(T) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texCSRJacobiDXref), LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
        }
    }
}

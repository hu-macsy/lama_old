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
 * @since 1.0.0
 */

#include <cuda_runtime.h>

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
    cudaStream_t stream );

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
    T* y_d );

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
    T* const y_d );

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
    T* y_d );

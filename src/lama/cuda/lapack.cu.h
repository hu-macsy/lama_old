/**
 * @file lapack.cu.h
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
 * @brief lapack.cu.h
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#include <cuda_runtime.h>

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
    cudaStream_t stream );

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
    const T omega );

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
    cudaStream_t stream );

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
    const T omega );

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
    cudaStream_t stream );

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
    const T omega );

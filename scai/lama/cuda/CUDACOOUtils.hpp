/**
 * @file CUDACOOUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Implementation of COO utilities with CUDA
 * @author Thomas Brandes
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/** This class provides CUDA parallelized routines needed for COO format.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDACOOUtils
{
public:

    /** Implementation for COOKernelTrait::offsets2ia with CUDA on GPUs */

    static void offsets2ia(
        IndexType cooIA[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );

    /** Implementation for COOKernelTrait::setCSRData with CUDA on GPUs */

    template<typename COOValueType,typename CSRValueType>
    static void setCSRData(
        COOValueType cooValues[],
        const CSRValueType csrValues[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );

    /** Implementation for COOKernelTrait::ia2offsets with CUDA on GPUs */

    static void ia2offsets(
        IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals,
        const IndexType cooIA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::normalGEMV with CUDA on GPUs */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[] );

    /** Implementation for COOKernelTrait::normalGEVM with CUDA on GPUs */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[] );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

    /** Constructor for registration. */

    CUDACOOUtils();

    /** Destructor for unregistration. */

    ~CUDACOOUtils();

    /** Static variable for registration at static initialization. */

    static CUDACOOUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */

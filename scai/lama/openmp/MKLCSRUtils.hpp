/**
 * @file MKLCSRUtils.hpp
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
 * @brief Implementation of CSR utilities with MKL for Host
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries

#include <scai/logging.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/** This class provides routines on compressed sparse row data
 */

class COMMON_DLL_IMPORTEXPORT MKLCSRUtils
{
public:

    /** Implementation for CSRKernelTrait::Transpose::convertCSR2CSC using MKL */

    template<typename ValueType>
    static void convertCSR2CSC(
        IndexType cscIA[],
        IndexType cscJA[],
        ValueType cscValues[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        IndexType numRows,
        IndexType numColumns,
        IndexType numValues );

    /** Implementation for CSRKernelTrait::Mult::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType nnz,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

    /** Routine that registers all routines of this class at the Kernel registry. */

    static void registerKernels();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    static bool initialized;

    static bool registerInterface();

};

} /* end namespace lama */

} /* end namespace scai */

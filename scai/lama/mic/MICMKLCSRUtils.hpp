/**
 * @file MICMKLCSRUtils.hpp
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
 * @brief Implementation of CSR utilities with MKL for MIC
 * @author Thomas Brandes
 * @date 02.07.2013
 * @since 1.1.0
 */
#ifndef LAMA_MIC_MKL_CSR_UTILS_HPP_
#define LAMA_MIC_MKL_CSR_UTILS_HPP_

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/LAMATypes.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

class SyncToken;
// forward declaration

/** This class provides routines on compressed sparse row data
 */

class COMMON_DLL_IMPORTEXPORT MICMKLCSRUtils
{
public:

    /** Implementation for CSRUtilsInterface::Mult::normalGEMV  */

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
        const ValueType csrValues[],
        SyncToken* syncToken );

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct CSRUtilsInterface& CSRUtils );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :

    static bool initialized;

    static bool registerInterface();

};

} /* end namespace lama */

} /* end namespace scai */

#endif //  LAMA_MIC_MKL_CSR_UTILS_HPP_

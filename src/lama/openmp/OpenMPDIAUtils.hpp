/**
 * @file OpenMPDIAUtils.hpp
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
 * @brief General conversion routines for DIA sparse matrices.
 * @author Thomas Brandes
 * @date 15.07.2012
 * $Id$
 */
#ifndef LAMA_OPENMP_DIA_UTILS_HPP_
#define LAMA_OPENMP_DIA_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

namespace lama
{

/** This class provides routines to converse DIA storage data to CSR storage data and vice versa.
 */

class LAMA_DLL_IMPORTEXPORT OpenMPDIAUtils
{
public:

    /** OpenMP implementation for DIAUtilsInterface::Counting::getCSRSizes */

    template<typename ValueType>
    static void getCSRSizes(
        IndexType csrSizes[],
        bool diagonalFlag,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[],
        const ValueType eps );

    /** OpenMP implementation for DIAUtilsInterface::Conversions::getCSRValues.  */

    template<typename DIAValueType,typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const bool diagonalFlag,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const DIAValueType diaValues[],
        const DIAValueType eps );

    /** Implementation for DIAUtilsInterface::Mult::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[],
        SyncToken* syncToken );

    /** Implementation for DIAUtilsInterface::Solver::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffset[],
        const ValueType diaValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega,
        const IndexType numRows,
        SyncToken* syncToken );

    /** Implemenatation for DIAUtilsInterface::Reductions::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[] );

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct DIAUtilsInterface& DIAUtils );

private:

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[] );

    static bool initialized;

    static bool registerInterface();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif //  LAMA_OPENMP_DIA_UTILS_HPP_

/**
 * @file OpenMPCOOUtils.hpp
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
 * @brief OpenMP implemenations for routines to be avaialble for COOUtilsInterface.
 * @author Thomas Brandes
 * @date 24.06.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/LAMATypes.hpp>
#include <scai/tasking/SyncToken.hpp>

// assert
#include <scai/lama/exception/LAMAAssert.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** This class provides OpenMP implementations to be used for COOUtilsInterface.
 *
 *  COOStorage is not well supported, but we provide conversions between COO
 *  and CSR as well as matrix times vector operation.
 */

class COMMON_DLL_IMPORTEXPORT OpenMPCOOUtils
{
public:

    /** OpenMP implementation for COOUtilsInterface::Counting::offsets2ia */

    static void offsets2ia(
        IndexType cooIA[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );

    /** OpenMP implementation for COOUtilsInterface::Counting::getCSRSizes */

    static void getCSRSizes(
        IndexType csrSizes[],
        const IndexType numRows,
        const IndexType numValues,
        const IndexType cooIA[] );

    /** OpenMP implementation for COOUtilsInterface::Conversions::getCSRValues */

    template<typename COOValueType,typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        IndexType csrIA[],
        const IndexType numRow,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const COOValueType cooValues[] );

    /** OpenMP implementation for COOUtilsInterface::Conversions::setCSRData */

    template<typename COOValueType,typename CSRValueType>
    static void setCSRData(
        COOValueType cooValues[],
        const CSRValueType csrValues[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );

    /** Implementation for CSRUtilsInterface::Mult::normalGEMV  */

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
        const ValueType cooValues[],
        tasking::SyncToken* syncToken );

    /** Implementation for CSRUtilsInterface::Mult::normalGEVM  */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        tasking::SyncToken* syncToken );

    /** Implementation for COOUtilsInterface::Solver::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType* solution,
        const IndexType cooNumValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega,
        const IndexType numRows,
        tasking::SyncToken* syncToken );

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct COOUtilsInterface& COOUtils );

private:

    static bool initialized;

    static bool registerInterface();

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */

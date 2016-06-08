/**
 * @file OpenMPCOOUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief OpenMP implemenations for routines to be avaialble for COOKernelTrait.
 * @author Thomas Brandes
 * @date 24.06.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides OpenMP implementations to be used for COOKernelTrait.
 *
 *  COOStorage is not well supported, but we provide conversions between COO
 *  and CSR as well as matrix times vector operation.
 */

class COMMON_DLL_IMPORTEXPORT OpenMPCOOUtils
{
public:

    /** OpenMP implementation for COOKernelTrait::hasDiagonalProperty */

    static bool hasDiagonalProperty (
        const IndexType cooIA[],
        const IndexType cooJA[],
        const IndexType n );

    /** OpenMP implementation for COOKernelTrait::offsets2ia */

    static void offsets2ia(
        IndexType cooIA[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );


    template<typename COOValueType,typename OtherValueType>
    static void scaleRows(
        COOValueType cooValues[],
        const OtherValueType rowValues[],
        const IndexType cooIA[],
        const IndexType numValues );
    
    /** OpenMP implementation for COOKernelTrait::getCSRSizes */

    static void getCSRSizes(
        IndexType csrSizes[],
        const IndexType numRows,
        const IndexType numValues,
        const IndexType cooIA[] );

    /** Serial implementation for COOKernelTrait::getCSRValues */

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

    /** OpenMP implementation for COOKernelTrait::setCSRData */

    template<typename COOValueType,typename CSRValueType>
    static void setCSRData(
        COOValueType cooValues[],
        const CSRValueType csrValues[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numDiagonals );

    /** Implementation for COOKernelTrait::normalGEMV  */

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

    /** Implementation for COOKernelTrait::normalGEVM  */

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
        const ValueType cooValues[] );

    /** Implementation for COOKernelTrait::jacobi  */

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
        const IndexType numRows );

private:

    template<typename ValueType>
    static void normalGEMV_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,
        const std::pair<ValueType, const ValueType*> by,
        const IndexType numRows,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[] );

    template<typename ValueType>
    static void normalGEVM_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,
        const std::pair<ValueType, const ValueType*> by,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[] );

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

    /** Constructor for registration. */

    OpenMPCOOUtils();

    /** Destructor for unregistration. */

    ~OpenMPCOOUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPCOOUtils guard;

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace sparsekernel */

} /* end namespace scai */

/**
 * @file OpenMPDIAUtils.hpp
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
 * @endlicense
 *
 * @brief General conversion routines for DIA sparse matrices.
 * @author Thomas Brandes
 * @date 15.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>


namespace scai
{

namespace sparsekernel
{

/** This class provides routines to converse DIA storage data to CSR storage data and vice versa.
 */

class COMMON_DLL_IMPORTEXPORT OpenMPDIAUtils
{
public:

    /** OpenMP implementation for DIAKernelTrait::getCSRSizes */

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

    /** OpenMP implementation for DIAKernelTrait::getCSRValues.  */

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

    /** Implementation for DIAKernelTrait::normalGEMV  */

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

    /** Implementation for DIAKernelTrait::normalGEVM  */

    template<typename ValueType>
    static void normalGEVM(
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

    /** Implementation for DIAKernelTrait::jacobi  */

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
        const IndexType numRows );

    /** Implemenatation for DIAKernelTrait::Reductions::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[] );

private:

    /** Help routine with max 9 arguments required */

    template<typename ValueType>
    static void normalGEMV_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,  // alpha, x
        const std::pair<ValueType, const ValueType*> by,  // beta, y
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[] );

    template<typename ValueType>
    static void normalGEVM_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,  // alpha, x
        const std::pair<ValueType, const ValueType*> by,  // beta, y
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[] );

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

    /** Constructor for registration. */

    OpenMPDIAUtils();

    /** Destructor for unregistration. */

    ~OpenMPDIAUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPDIAUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace sparsekernel */

} /* end namespace scai */

/**
 * @file OpenMPCOOUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
#include <scai/common/MatrixOp.hpp>

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

    /** OpenMP implementation for COOKernelTrait::getValuePos */

    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const IndexType numValues );

    /** OpenMP implementation for COOKernelTrait::offsets2ia */

    static void offsets2ia(
        IndexType cooIA[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows );

    /** OpenMP implementation for COOKernelTrait::ia2offsets */

    static void ia2offsets(
        IndexType csrIA[],
        const IndexType numRows,
        const IndexType cooIA[],
        const IndexType numValues );

    template<typename ValueType>
    static void scaleRows(
        ValueType cooValues[],
        const ValueType rowValues[],
        const IndexType cooIA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const common::MatrixOp op );

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

    /** Implementation for COOKernelTrait::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType cooNumValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const ValueType localDiagonal[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numRows );

    /** Implementation for COOKernelTrait::hasDiagonalProperty */

    static bool hasDiagonalProperty(
        const IndexType numDiagonals,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::getDiagonal */

    template<typename ValueType>
    static void getDiagonal(
        ValueType diagonal[],
        const IndexType numDiagonals,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::setDiagonalV */

    template<typename ValueType>
    static void setDiagonalV(
        ValueType cooValues[],
        const ValueType diagonal[],
        const IndexType numDiagonals,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::setDiagonal */

    template<typename ValueType>
    static void setDiagonal(
        ValueType cooValues[],
        const ValueType diagonal,
        const IndexType numDiagonals,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::getRow */

    static IndexType getRow(
        IndexType& offset,
        const IndexType cooIA[],
        const IndexType numValues,
        const IndexType i );

    /** Implementation for COOKernelTrait::getColumn */

    static IndexType getColumn(
        IndexType positions[],
        const IndexType cooJA[],
        const IndexType numValues,
        const IndexType j );

private:

    static IndexType getRowStartPos( const IndexType i,
                                     const IndexType cooIA[],
                                     const IndexType numValues );

    /** Struct for registration of methods without template arguments */

    struct Registrator
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

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

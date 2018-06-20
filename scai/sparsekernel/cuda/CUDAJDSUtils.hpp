/**
 * @file CUDAJDSUtils.hpp
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
 * @brief Implementation of JDS utilities with CUDA
 * @author Thomas Brandes
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>
#include <scai/common/BinaryOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides CUDA parallelized routines needed for JDS format.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDAJDSUtils
{
public:

    /** CUDA implementation of JDSKernelTrait::setRows */

    template<typename ValueType>
    static void setRows(
        ValueType jdsValues[],
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        const ValueType rowValues[],
        common::BinaryOp op );

    /** CUDA implementation of JDSKernelTrait::getRow */

    template<typename ValueType>
    static void getRow(
        ValueType row[],
        const IndexType i,
        const IndexType numColumns,
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType ja[],
        const ValueType values[] );

    /** CUDA implementation of JDSKernelTrait::getValuePos */

    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType perm[],
        const IndexType ja[] );

    /** CUDA implementation for JDSKernelTrait::ilg2dlg */

    static IndexType ilg2dlg(
        IndexType dlg[],
        const IndexType numDiagonals,
        const IndexType ilg[],
        const IndexType numRows );

    /** Conversion of JDS to CSR as specified in JDSKernelTrait::getCSRValues  */

    template<typename JDSValueType, typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType jdsPerm[],
        const IndexType jdsILG[],
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const JDSValueType jdsValues[] );

    /** Conversion of CSR to JDS in CUDA as specified in JDSKernelTrait::getCSRValues  */

    template<typename JDSValueType, typename CSRValueType>
    static void setCSRValues(
        IndexType jdsJA[],
        JDSValueType jdsValues[],
        const IndexType numRows,
        const IndexType jdsPerm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** Implementation for JDSKernelTrait::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType numRows,
        const IndexType jdsPerm[],
        const IndexType jdsIlg[],
        const IndexType jdsNumDiagonals,
        const IndexType jdsDlg[],
        const IndexType jdsJA[],
        const ValueType jdsValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega );

    /** Implementation for JDSKernelTrait::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solutionLocal[],
        const IndexType numRows,
        const ValueType diagonal[],
        const IndexType ndlg_halo,
        const IndexType jdsPermHalo[],
        const IndexType jdsIlgHalo[],
        const IndexType jdsDlgHalo[],
        const IndexType jdsJAHalo[],
        const ValueType jdsValuesHalo[],
        const ValueType oldSolutionHalo[],
        const ValueType omega );

    /** CUDA implementation for JDSKernelTrait::normalGEMV */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType perm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[],
        const common::MatrixOp op );

    /** Implementation for JDSKernelTrait::sparseGEMV with CUDA on GPU */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType perm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[],
        const common::MatrixOp op );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    /** Struct for registration of methods with two template arguments.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for all supported types.
     */

    template<typename ValueType, typename OtherValueType>
    struct RegistratorVO
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDAJDSUtils();

    /** Destructor for unregistration. */

    ~CUDAJDSUtils();

    /** Static variable for registration at static initialization. */

    static CUDAJDSUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */

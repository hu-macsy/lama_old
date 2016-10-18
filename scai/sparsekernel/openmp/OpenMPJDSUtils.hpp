/**
 * @file OpenMPJDSUtils.hpp
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
 * @brief Implementation of JDS utilities with OpenMP
 * @author Thomas Brandes
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/BinaryOp.hpp>

#include <utility>

namespace scai
{

namespace sparsekernel
{

/** This class provides OpenMP implementations as needed for JDSUtilKernelTrait.  */

class COMMON_DLL_IMPORTEXPORT OpenMPJDSUtils
{
public:

    /** This method scales the matrix using an value vector */

    template<typename ValueType, typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        ValueType mValues[],
        const OtherValueType values[] );

    /** This method sets row as dens vector of the i'th row of the matrix */

    template<typename ValueType, typename OtherValueType>
    static void getRow(
        OtherValueType row[],
        const IndexType i,
        const IndexType numColumns,
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType ja[],
        const ValueType values[] );

    template<typename ValueType, typename OtherValueType>
    static void setRow(
        ValueType values[],
        const IndexType i,
        const IndexType numColumns,
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType ja[],
        const OtherValueType row[],
        const utilskernel::binary::BinaryOp op );

    /** Implementation for JDSKernelTrait::getValuePos */

    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType perm[],
        const IndexType ja[] );

    /** Implementation for JDSKernelTrait::getValuePosCol */

    static IndexType getValuePosCol(
        IndexType row[],
        IndexType pos[],
        const IndexType j,
        const IndexType numRows,
        const IndexType ilg[],
        const IndexType dlg[],
        const IndexType perm[],
        const IndexType ja[] );

    /** This method checks if the matrix has diagonal property */

    static bool checkDiagonalProperty(
        const IndexType numDiagonals,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType perm[],
        const IndexType ja[],
        const IndexType dlg[] );

    /** Bucket sort for row sorting */

    static void sortRows( IndexType array[], IndexType perm[], const IndexType n );

    /** Compute the inverse permutation as specified in JDSUtilKernelTrait::Sort::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** Compute dlg array from ilg array as specified in JDSUtilKernelTrait::Conversions::ilg2dlg */

    static IndexType ilg2dlg(
        IndexType dlg[],
        const IndexType numDiagonals,
        const IndexType ilg[],
        const IndexType numRows );

    /** Conversion of JDS to CSR as specified in JDSUtilKernelTrait::Conversions::getCSRValues  */

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

    /** Conversion of CSR to JDS as specified in JDSUtilKernelTrait::Conversions::setCSRValues. */

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

    /** Implementation for JDSUtilKernelTrait::Mult:normalGEMV with OpenMP on Host */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType perm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[] );

    /** Implementation for JDSUtilKernelTrait::Mult:normalGEVM with OpenMP on Host */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numColumns,
        const IndexType perm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[] );

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType numRows,
        const IndexType jdsPerm[],
        const IndexType jdsILG[],
        const IndexType jdsNumDiagonals,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega );

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType numRows,
        const ValueType localDiagonal[],
        const IndexType numDiagonals,
        const IndexType jdsHaloPerm[],
        const IndexType jdsHaloILG[],
        const IndexType jdsHaloDLG[],
        const IndexType jdsHaloJA[],
        const ValueType jdsHaloValues[],
        const ValueType oldSolution[],
        const ValueType omega );

private:

    // We need for asynchronous execution versions with max 9 args

    template<typename ValueType>
    static void normalGEMV_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,     // alpha, x
        const std::pair<ValueType, const ValueType*> by,     // beta, y
        const std::pair<IndexType, const IndexType*> rows,   // nrows, jdsILG,
        const IndexType perm[],
        const std::pair<IndexType, const IndexType*> dlg,    // ndlg, jdsDLG
        const IndexType jdsJA[],
        const ValueType jdsValues[] );

    template<typename ValueType>
    static void normalGEVM_a(
        ValueType result[],
        const std::pair<ValueType, const ValueType*> ax,     // alpha, x
        const std::pair<ValueType, const ValueType*> by,     // beta, y
        const std::pair<IndexType, const IndexType*> rows,   // nrows, jdsILG,
        const IndexType perm[],
        const std::pair<IndexType, const IndexType*> dlg,    // ndlg, jdsDLG
        const IndexType jdsJA[],
        const ValueType jdsValues[] );

    /** Routine that registers all methods at the kernel registry. */

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

    OpenMPJDSUtils();

    /** Destructor for unregistration. */

    ~OpenMPJDSUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPJDSUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */

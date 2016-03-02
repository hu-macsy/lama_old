/**
 * @file OpenMPJDSUtils.hpp
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
 * @brief Implementation of JDS utilities with OpenMP
 * @author Thomas Brandes
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/Scalar.hpp>

// internal scai libraries
#include <scai/kregistry/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

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

    template<typename ValueType,typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType perm[],
        const IndexType ilg[],
        const IndexType dlg[],
        ValueType mValues[],
        const OtherValueType values[] );

    /** This method sets row as dens vector of the i'th row of the matrix */

    template<typename ValueType,typename OtherValueType>
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

    template<typename ValueType>
    static ValueType getValue(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType* dlg,
        const IndexType* ilg,
        const IndexType* perm,
        const IndexType* ja,
        const ValueType* values );

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

    template<typename JDSValueType,typename CSRValueType>
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

    template<typename JDSValueType,typename CSRValueType>
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

    /** Routine that registers all methods at the kernel registry. */

    SCAI_DECLARE_REGISTRATOR( Registrator )
    SCAI_DECLARE_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_DECLARE_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

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

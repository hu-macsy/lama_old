/**
 * @file CUDAJDSUtils.hpp
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
 * @brief Implementation of JDS utilities with CUDA
 * @author Thomas Brandes
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
//#include <scai/lama/Scalar.hpp>

// internal scai library
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

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

    /** Implementation for function type JDSKernelTrait::sortRows  */

    static void sortRows( IndexType array[], IndexType perm[], const IndexType n );

    /** Compute dlg array from ilg array as specified in JDSKernelTrait::ilg2dlg */

    static IndexType ilg2dlg(
        IndexType dlg[],
        const IndexType numDiagonals,
        const IndexType ilg[],
        const IndexType numRows );

    /** Compute the inverse permutation as specified in JDSKernelTrait::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** Conversion of JDS to CSR as specified in JDSKernelTrait::getCSRValues  */

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

    /** Conversion of CSR to JDS in CUDA as specified in JDSKernelTrait::getCSRValues  */

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

    /** Implementation for JDSKernelTrait::normalGEMV with CUDA on GPU */

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

    /** Implementation for JDSKernelTrait::normalGEVM with CUDA on GPU */

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
        const ValueType jdsValues[] );

    /** Implementation for JDSKernelTrait::sparseGEVM with CUDA on GPU */

    template<typename ValueType>
    static void sparseGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType perm[],
        const IndexType jdsILG[],
        const IndexType ndlg,
        const IndexType jdsDLG[],
        const IndexType jdsJA[],
        const ValueType jdsValues[] );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

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

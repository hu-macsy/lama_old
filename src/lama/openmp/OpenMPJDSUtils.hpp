/**
 * @file OpenMPJDSUtils.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */
#ifndef LAMA_OPENMP_JDS_UTILS_HPP_
#define LAMA_OPENMP_JDS_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>
#include <lama/Scalar.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** This class provides OpenMP implementations as needed for JDSUtilsInterface.  */

class LAMA_DLL_IMPORTEXPORT OpenMPJDSUtils
{
public:

    /** This method replaces the diagonal */

    template<typename ValueType>
    static void setDiagonalWithScalar( const IndexType numDiagonal, ValueType values[], Scalar scalar );

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

    template<typename ValueType,typename NoType>
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

    /** This method checks if the matrix indexes are correct */

    static bool check(
        const IndexType numRows,
        const IndexType numValues,
        const IndexType numColumns,
        const IndexType ja[],
        const IndexType ilg[],
        const IndexType dlg[] );

    /** Bucket sort for row sorting */

    static void sortRows( IndexType array[], IndexType perm[], const IndexType n );

    /** Compute the inverse permutation as specified in JDSUtilsInterface::Sort::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** Compute dlg array from ilg array as specified in JDSUtilsInterface::Conversions::ilg2dlg */

    static IndexType ilg2dlg(
        IndexType dlg[],
        const IndexType numDiagonals,
        const IndexType ilg[],
        const IndexType numRows );

    /** Conversion of JDS to CSR as specified in JDSUtilsInterface::Conversions::getCSRValues  */

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

    /** Conversion of CSR to JDS as specified in JDSUtilsInterface::Conversions::setCSRValues. */

    template<typename JDSValueType,typename CSRValueType>
    static void setCSRValues(
        IndexType jdsJA[],
        JDSValueType jdsValues[],
        const IndexType numRows,
        const IndexType jdsPerm[],
        const IndexType jdsILG[],
        const IndexType jdsDLG[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** Implementation for JDSUtilsInterface::Mult:normalGEMV with OpenMP on Host */

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
        const ValueType jdsValues[],
        SyncToken* syncToken );

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
        const ValueType omega,
        SyncToken* syncToken );

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType numRows,
        const ValueType localDiagonal[],
        const IndexType numDiagonals,
        const IndexType jdsHaloPerm[],
        const IndexType jdsHaloIlg[],
        const IndexType jdsHaloDlg[],
        const IndexType jdsHaloJA[],
        const ValueType jdsHaloValues[],
        const ValueType oldSolution[],
        const ValueType omega,
        SyncToken* syncToken );

    /** Method for registration of module routines at the interface. */

    static void setInterface( struct JDSUtilsInterface& JDSUtils );

private:

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

};

/* --------------------------------------------------------------------------- */

} // namespace lama

#endif //  LAMA_OPENMP_JDS_UTILS_HPP_

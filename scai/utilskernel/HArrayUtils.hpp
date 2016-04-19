/*
 * @file HArrayUtils.hpp
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
 * @brief Definition of class with utility routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>
#include <scai/common/ReductionOp.hpp>
#include <scai/common/mepr/TemplateSpecifier.hpp>

namespace scai
{

namespace utilskernel
{

/** Class that contains some utility routines used at several places. */

class COMMON_DLL_IMPORTEXPORT HArrayUtils
{

public:

    /** Static method for HArray assignment with type conversions.
     *
     *  @param[out] target   contains copy of source values
     *  @param[in]  source   array with source values
     *  @param[in]  prefLoc  specifies optionally at which preferred context target will have valid values
     *
     *  \code
     *  HArray<float> fvalues;
     *  HArray<double> dvalues;
     *  ...
     *  dvalues = fvalues;   // not supported
     *  assign( dvalues, fvalues );  // supported
     *  \endcode
     *  Size of target array will be the same as the source array.
     */
    static void assign( 
        hmemo::_HArray& target, 
        const hmemo::_HArray& source, 
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    static void assignOp( 
        hmemo::_HArray& target, 
        const hmemo::_HArray& source, 
        const common::reduction::ReductionOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** 
     *  @brief Gathering (unstructured read) of values with heterogeneous arrays. 
     *
     *  target[i] = source[index[i]] 
     */
    static void assignGather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& index,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read of values) with HArrays, template typed version
     */
    template<typename TargetValueType,typename SourceValueType>
    static void gather(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::HArray<IndexType>& index,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** 
     *  @brief Scattering (unstructured write) of values with heterogeneous arrays. 
     *
     *  target[index[i]] = source[i] 
     */
    static void assignScatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& index,
        const hmemo::_HArray& source,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Scatter (unstructured write of values) with HArrays, template typed version
     */
    template<typename TargetValueType,typename SourceValueType>
    static void scatter(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<IndexType>& index,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Setting one scalar element for all elements of HArray
     *
     *  target[i] <op>= value 
     */
    template<typename ValueType>
    static void assignScalar( 
        hmemo::_HArray& target,
        const ValueType value,
        const common::reduction::ReductionOp op, 
        hmemo::ContextPtr prefLoc  = hmemo::ContextPtr() )
        __attribute__( ( noinline ) );

    /** This method sets a single value in a heterogeneous array.
     *
     *  @param[in,out] target  Harray where a value to set
     *  @param[in]     index  position to set ( 0 <= index < target.size() )
     *  @param[in]     val    value to set
     *
     *  The value will be set at a valid context.
     */

    template<typename ValueType>
    static void setVal( hmemo::_HArray& target, const IndexType index, const ValueType val );

    template<typename ValueType>
    static ValueType getVal( const hmemo::_HArray& array, const IndexType index );

    /** Scaled assignment: result = beta * y
     *
     *  @param[out] result  output array
     *  @param[in]  beta    scaling factor
     *  @param[in]  y       source array
     *  @param[in]  context location where operation is done
     */

    template<typename ValueType>
    static void assignScaled(
        hmemo::HArray<ValueType>& result,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Axpy: result += beta * y
     *
     *  @param[in,out] result  output array
     *  @param[in]     beta    scaling factor
     *  @param[in]     y       source array
     *  @param[in]     prefLoc location where operation is done
     */

    template<typename ValueType>
    static void axpy(
        hmemo::HArray<ValueType>& result,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Addition of two arrays: result = alpha * x + beta * y
     *
     *  @param[out] result  output array
     *  @param[in]  alpha   scaling factor
     *  @param[in]  x       source array
     *  @param[in]  beta    scaling factor
     *  @param[in]  y       source array
     *  @param[in]  context location where operation is done
     */

    template<typename ValueType>
    static void arrayPlusArray(    
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** scale array in place : array *= beta
     *
     *  Note: scale will be done where array has currently valid values. The preferred
     *        location is not taken if the array is not valid there.
     */

    template<typename ValueType>
    static void scale( hmemo::HArray<ValueType>& array, const ValueType beta, hmemo::ContextPtr prefLoc );

    /** Replace in a complex array its values with the conjugate values */

    template<typename ValueType>
    static void conj( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /*
     * Implementation of functions
     */
    template<typename TargetValueType,typename SourceValueType>
    static void setArray(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const common::reduction::ReductionOp op,
        hmemo::ContextPtr context );

    template<typename ValueType>
    static void setScalar(
        hmemo::HArray<ValueType>& target,
        const ValueType value,
        const common::reduction::ReductionOp op, 
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() )
        __attribute__( ( noinline ) );

    template<typename ValueType>
    static void setValImpl(
        hmemo::HArray<ValueType>& target,
        const IndexType index,
        const ValueType val );

    template<typename ValueType>
    static ValueType getValImpl(
        const hmemo::HArray<ValueType>& array,
        const IndexType index );

    template<typename ValueType>
    static ValueType reduce( 
        const hmemo::HArray<ValueType>& array,
        const common::reduction::ReductionOp redOp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType asum( 
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType absMaxDiffVal( 
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType dotProduct( 
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** array = 1.0 / array elementwise */

    template<typename ValueType>
    static void invert( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check for an index array whether all values are smaller than n */

    static bool validIndexes( const hmemo::HArray<IndexType>& array, const IndexType size, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check whether values in array are sorted ascending or descending. */

    template<typename ValueType>
    static bool isSorted( const hmemo::HArray<ValueType>& array, const bool isAscending, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build the running sums for an array; note that result array will contain one element more. 
     *
     *  @param[in,out]  array contains values for which running sum is built
     *  @param[in]      prefLoc optional the context where computation should be done
     *
     *  \code
     *       array( in ) = { 3, 5, 7, 2 }, array( out ) = { 0, 3, 8, 15, 17 }
     *  \endcode
     */

    template<typename ValueType>
    static ValueType scan( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static void sort( hmemo::HArray<ValueType>& array, hmemo::HArray<IndexType>& perm, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Initialize an array with the sequence 0, .., n-1 
     *
     *  @param[out] array will contain the values 0, ..., n-1
     *  @param[in]  n     becomes size of the array
     *  @param[in]  prefLoc optional the context where allocation/initialization should be done
     */

    static void setOrder( hmemo::HArray<IndexType>& array, IndexType n, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Sete an array with random values. 
     *
     *  @param[out] array    will contain random values of its type
     *  @param[in]  n        number of values, becomes size of array
     *  @param[in]  fillRate ratio of non-zero values
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */

    template<typename ValueType>
    static void setRandom( hmemo::HArray<ValueType>& array, 
                           IndexType n, 
                           float fillRate = 1.0f, 
                           hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build sparse array from dense array, needed for conversion DenseVector -> SparseVector */

    template<typename ValueType>
    static void buildSparseArray( 
        hmemo::HArray<ValueType>& sparseArray,
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::HArray<ValueType>& denseArray,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build dense array from sparse array, needed for conversion SparseVector -> DenseVector 
     *
     *  @param[out]  denseArray will contain the dense data (input its only its size)
     *  @param[in]   denseN is the size of the dense array
     *  @param[in]   sparseArray contains non-zero values
     *  @param[in]   sparseIndexes are the positions of the non-zero values
     *  @param[in]   prefLoc is the context where operation should be done
     *
     *  Note: sparseIndexes must contain only indexes between 0 and denseN - 1
     */

    template<typename ValueType>
    static void buildDenseArray( 
        hmemo::HArray<ValueType>& denseArray,
        const IndexType denseN,
        const hmemo::HArray<ValueType>& sparseArray,
        const hmemo::HArray<IndexType>& sparseIndexes,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    SCAI_DECLARE_TEMPLATESPECIFIER( SpecifierV, template<typename ValueType> )

    HArrayUtils();  // static class, no objects outside

    static HArrayUtils guard;   // dummy object guarantees template method instantiation
};

} /* end namespace utilskernel */

} /* end namespace scai */

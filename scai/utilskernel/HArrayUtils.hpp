/**
 * @file HArrayUtils.hpp
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
 * @brief Definition of class with utility routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>
#include <scai/utilskernel/ReductionOp.hpp>
#include <scai/utilskernel/ElementwiseOp.hpp>

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
        const reduction::ReductionOp op,
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
    template<typename TargetValueType, typename SourceValueType>
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
        const reduction::ReductionOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Scatter (unstructured write of values) with HArrays, template typed version
     * 
     *  target[index[i]] = source[i]
     */
    template<typename TargetValueType, typename SourceValueType>
    static void scatter(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<IndexType>& index,
        const hmemo::HArray<SourceValueType>& source,
        const reduction::ReductionOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Setting one scalar element for all elements of HArray
     *
     *  target[i] _op_= value
     */
    template<typename ValueType>
    static void assignScalar(
        hmemo::_HArray& target,
        const ValueType value,
        const reduction::ReductionOp op,
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
     *  @param[in]  prefLoc location where operation should be done if possible
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
     *  @param[in]  prefLoc location where operation should be done if possible
     */

    template<typename ValueType>
    static void arrayPlusArray(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Multiplication of two arrays: result = alpha * x * y
     *
     *  @param[out] result  output array
     *  @param[in]  alpha   scalar that is multiplied for each element
     *  @param[in]  x       first input array
     *  @param[in]  y       second input array
     *  @param[in]  prefLoc location where operation should be done if possible
     *
     *  Any alias of the arrays x, y, result is supported.
     */

    template<typename ValueType>
    static void arrayTimesArray(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** scale array in place : array *= beta
     *
     *  Note: scale will be done where array has currently valid values. The preferred
     *        location is not taken if the array is not valid there.
     */

    template<typename ValueType>
    static void scale( hmemo::HArray<ValueType>& array, const ValueType beta, hmemo::ContextPtr prefLoc );

    /*
     * Implementation of functions
     */
    template<typename TargetValueType, typename SourceValueType>
    static void setArray(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const reduction::ReductionOp op,
        hmemo::ContextPtr context );

    /** General version for setting sectioned arrays. 
     *
     *  @param target is the target array
     *  @param targetOffset is the offset in the target array
     *  @param targetStride is the stride used in target array
     *  @param source is the source array
     *  @param sourceOffset is the offset in the source array
     *  @param sourceStride is the stride used in source array
     *  @param n is the number of elements to set
     *  @param op specifies how to combine old and new value
     */

    template<typename TargetValueType, typename SourceValueType>
    static void setArraySection(
        hmemo::HArray<TargetValueType>& target,
        const IndexType targetOffset,
        const IndexType targetStride,
        const hmemo::HArray<SourceValueType>& source,
        const IndexType sourceOffset,
        const IndexType sourceStride,
        const IndexType n,
        const reduction::ReductionOp op,
        hmemo::ContextPtr context );

    template<typename ValueType>
    static void setScalar(
        hmemo::HArray<ValueType>& target,
        const ValueType value,
        const reduction::ReductionOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() )
    __attribute__( ( noinline ) );

    template<typename ValueType>
    static void setValImpl(
        hmemo::HArray<ValueType>& target,
        const IndexType index,
        const ValueType val,
        const reduction::ReductionOp op );

    template<typename ValueType>
    static ValueType getValImpl(
        const hmemo::HArray<ValueType>& array,
        const IndexType index );

    template<typename ValueType>
    static ValueType reduce(
        const hmemo::HArray<ValueType>& array,
        const reduction::ReductionOp redOp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType asum(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType nrm2(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType absMaxDiffVal(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static void copysign(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType dotProduct(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** execute elementwise functions */

    template<typename ValueType>
    static void execElementwise(
        hmemo::HArray<ValueType>& array,
        const elementwise::ElementwiseOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** execute elementwise pow with base and exponent array */

    template<typename ValueType>
    static void pow( 
        hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** execute elementwise pow with base array and */

    template<typename ValueType>
    static void powBase( 
        hmemo::HArray<ValueType>& array,
        const ValueType base,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static void powExp( 
        hmemo::HArray<ValueType>& array,
        const ValueType exp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check for an index array whether all values are smaller than n */

    static bool validIndexes( const hmemo::HArray<IndexType>& array, const IndexType size, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check whether values in array are sorted ascending or descending. */

    template<typename ValueType>
    static bool isSorted( const hmemo::HArray<ValueType>& array, const bool isAscending, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build the running sums for an array; note that result array will contain one element more.
     *
     *  @param[in,out]  array contains values for which running sum is built
     *  @param[in]      prefLoc optional the context where computation should be done
     *  @returns        the total sum and last value in the array is returned
     *
     *  \code
     *       array( in ) = { 3, 5, 7, 2 }, array( out ) = { 0, 3, 8, 15, 17 }
     *  \endcode
     */

    template<typename ValueType>
    static ValueType scan( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build the differences, opposite to scan, especially used to convert offsets to sizes
     *
     *  @param[in,out]  array contains values for which differences are built
     *  @param[in]      prefLoc optional the context where computation should be done
     *  @returns        value of the first element
     *
     *  \code
     *       array( in ) = { 0, 3, 8, 15, 17 }, array( out ) = { 3, 5, 7, 2 }, returns 0
     *       array( in ) = { 1, 3, 8, 15, 17 }, array( out ) = { 2, 5, 7, 2 }, returns 1
     *  \endcode
     */

    template<typename ValueType>
    static ValueType unscan( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Sort an array of values
     *
     *  @param[in,out] array is the array of values to be sorted
     *  @param[out] perm is the permutation that gives the sorted array
     *  @param[in] prefLoc is the preferred context where computation should be done
     * 
     *  Note: array_out = array_in[ perm ]
     *
     *  ToDo: ascending or descending, why no choice
     */

    template<typename ValueType>
    static void sort( hmemo::HArray<ValueType>& array, hmemo::HArray<IndexType>& perm, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Bucket sort of an array with integer values
     *
     *  @param[out] perm is permutation to get array sorted
     *  @param[out] offsets is an offset array for perm to sort it bucketwise, size is nb + 1
     *  @param[in] nb is the number of buckets
     *  @param[in] array contains bucket indexes, 0 <= array[i] < nb
     *  @param[in] prefLoc is the preferred context where computation should be done
     *
     *  Note: the sorted array is given by array[perm] 
     *  Note: perm.size() == array.size() if all values of array are correct bucket indexes
     *        otherwise perm.size() < array.size(), can still be used to get legal sorted buckets
     *  Note: in contrary to sort the array remains unchanged
     */

    template<typename BucketType>
    static void bucketSort(
        hmemo::HArray<IndexType>& offsets,
        hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<BucketType>& array,
        const BucketType nb,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Slighter version of bucketSort, counts only values */

    template<typename BucketType>
    static void bucketCount(
        hmemo::HArray<IndexType>& bucketSizes,
        const hmemo::HArray<BucketType>& array,
        const BucketType nb,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Initialize an array with the sequence 0, .., n-1
     *
     *  @param[out] array   will contain the values 0, ..., n-1
     *  @param[in]  n       becomes size of the array
     *  @param[in]  prefLoc optional the context where allocation/initialization should be done
     */

    static void setOrder( hmemo::HArray<IndexType>& array, IndexType n, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Initilize an array with a sequence of values starting with startValue, incrementing by inc
     *
     *  @param[out] array       will contain the values startValue, ..., startValue + (n-1)*inc
     *  @param[in]  startValue  startValue of the sequence
     *  @param[in]  inc         increment of the sequence
     *  @param[in]  n           becomes size of the array
     *  @param[in]  prefLoc     optional the context where allocation/initialization should be done
     */

    template<typename ValueType>
    static void setSequence( hmemo::HArray<ValueType>& array, ValueType startValue, ValueType inc, IndexType n, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Set an array with random values.
     *
     *  @param[out] array    will contain random values of its type
     *  @param[in]  n        number of values, becomes size of array
     *  @param[in]  fillRate ratio of non-zero values
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */

    static void setRandom( hmemo::_HArray& array,
                           IndexType n,
                           float fillRate = 1.0f,
                           hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Sete an array with random values.
     *
     *  @param[out] array    will contain random values of its type
     *  @param[in]  n        number of values, becomes size of array
     *  @param[in]  fillRate ratio of non-zero values
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */

    template<typename ValueType>
    static void setRandomImpl( hmemo::HArray<ValueType>& array,
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

    HArrayUtils();  // static class, no objects outside
    HArrayUtils( const HArrayUtils& );

};

} /* end namespace utilskernel */

} /* end namespace scai */

/**
 * @file HArrayUtils.hpp
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
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>

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

    /**
     *  @brief Gathering (unstructured read) of values with heterogeneous arrays.
     *
     *  target[i] op= source[indexes[i]]
     */
    static void gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    static void sparseGather(
        hmemo::_HArray& target,
        const hmemo::_HArray& sourceNonZeroValues,
        const hmemo::HArray<IndexType>& sourceNonZeroIndexes,
        const hmemo::HArray<IndexType>& indexes,
        const common::binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read of values) with HArrays, template typed version
     */
    template<typename TargetValueType, typename SourceValueType>
    static void gatherImpl(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read of values) for sparse data, typed version
     */
    template<typename TargetValueType, typename SourceValueType>
    static void sparseGatherImpl(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& sourceNonZeroValues,
        const hmemo::HArray<IndexType>& sourceNonZeroIndexes,
        const hmemo::HArray<IndexType>& indexes,
        const common::binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Scattering (unstructured write) of values with heterogeneous arrays.
     *
     *  target[index[i]] = source[i]
     */
    static void scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& index,
        const bool unique,
        const hmemo::_HArray& source,
        const common::binary::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Scatter (unstructured write of values) with HArrays, template typed version
     *
     *  target[index[i]] = source[i]
     */
    template<typename TargetValueType, typename SourceValueType>
    static void scatterImpl(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<IndexType>& index,
        const bool unique,
        const hmemo::HArray<SourceValueType>& source,
        const common::binary::BinaryOp op,
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
        const common::binary::BinaryOp op,
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

    /** Addition of scaled array with scalar: result = alpha * x + beta (elementwise)
     *
     *  @param[out] result  output array
     *  @param[in]  alpha   scaling factor
     *  @param[in]  x       source array
     *  @param[in]  beta    scaling factor
     *  @param[in]  prefLoc location where operation should be done if possible
     */

    template<typename ValueType>
    static void arrayPlusScalar(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Element-wise multiplication of two arrays: result[i] = alpha * x[i] * y[i]
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

    /*
     * Implementation of functions
     */
    static void setArray(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    template<typename TargetValueType, typename SourceValueType>
    static void setArrayImpl(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

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
     *  @param context is the preferred context where option is done
     */

    static void setArraySection(
        hmemo::_HArray& target,
        const IndexType targetOffset,
        const IndexType targetStride,
        const hmemo::_HArray& source,
        const IndexType sourceOffset,
        const IndexType sourceStride,
        const IndexType n,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    /** Typed version for setting sectioned arrays */

    template<typename TargetValueType, typename SourceValueType>
    static void setArraySectionImpl(
        hmemo::HArray<TargetValueType>& target,
        const IndexType targetOffset,
        const IndexType targetStride,
        const hmemo::HArray<SourceValueType>& source,
        const IndexType sourceOffset,
        const IndexType sourceStride,
        const IndexType n,
        const common::binary::BinaryOp op = common::binary::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    template<typename ValueType>
    static void setScalar(
        hmemo::HArray<ValueType>& target,
        const ValueType value,
        const common::binary::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() )
    __attribute__( ( noinline ) );

    template<typename ValueType>
    static void setValImpl(
        hmemo::HArray<ValueType>& target,
        const IndexType index,
        const ValueType val,
        const common::binary::BinaryOp op );

    template<typename ValueType>
    static ValueType getValImpl(
        const hmemo::HArray<ValueType>& array,
        const IndexType index );

    template<typename ValueType>
    static ValueType reduce(
        const hmemo::HArray<ValueType>& array,
        const common::binary::BinaryOp redOp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType reduce2(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        const common::binary::BinaryOp binOp,
        const common::binary::BinaryOp redOp,
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
    static ValueType dotProduct(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise unary operation on array: result[i] = op( x[i] )
     *
     *  @param[out] result  output array
     *  @param[in]  x       input array
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     */

    template<typename ValueType>
    static void unaryOp(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const common::unary::UnaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise binary operation on array: result[i] = op( x[i], y[i] )
     *
     *  @param[out] result  output array
     *  @param[in]  x       input array
     *  @param[in]  y       input array
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     */
    template<typename ValueType>
    static void binaryOp(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const hmemo::HArray<ValueType>& y,
        const common::binary::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise binary operation on array: result[i] = op( x, y[i] ), first arg is scalar
     *
     *  @param[out] result  output array
     *  @param[in]  x       input value
     *  @param[in]  y       input array
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     */
    template<typename ValueType>
    static void binaryOpScalar1(
        hmemo::HArray<ValueType>& result,
        const ValueType x,
        const hmemo::HArray<ValueType>& y,
        const common::binary::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise binary operation on array: result[i] = op( x[i], y ), second arg is scalar
     *
     *  @param[out] result  output array
     *  @param[in]  x       input array
     *  @param[in]  y       input value
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     *
     *  Note: this operation is different to binaryOpScalar1( result, y, x, loc ) if op is not commutative
     */
    template<typename ValueType>
    static void binaryOpScalar2(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const ValueType y,
        const common::binary::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check for an index array whether all values are smaller than n */

    static bool validIndexes( const hmemo::HArray<IndexType>& array, const IndexType size, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check whether values in array are sorted  
     *
     *  @param[in] array
     *  @param[in] op
     *  @param[in] prefloc 
     *  @returns   true if a[i] op a[i+1] is true for all neighbored pairs
     */
    template<typename ValueType>
    static bool isSorted( 
        const hmemo::HArray<ValueType>& array, 
        const common::binary::CompareOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build the running sums within an array
     *
     *  @param[in,out]  array contains values for which running sum is built
     *  @param[in]      first is an offset that is used as starting value
     *  @param[in]      exclusive if set the current values is not included in the running sum
     *  @param[in]      prefLoc optional the context where computation should be done
     *  @returns        the total sum and last value in the array is returned
     *
     *  \code
     *       array( in ) = { 3, 5, 7, 2 }, array( out ) = { 13, 18, 25, 27 }  first = 10, exclusive = false
     *       array( in ) = { 3, 5, 7, 2 }, array( out ) = { 10, 13, 18, 25 }  first = 10, exclusive = true
     *  \endcode
     */
    template<typename ValueType>
    static ValueType scan( 
        hmemo::HArray<ValueType>& array, 
        const ValueType first,
        bool exclusive, 
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

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
    static ValueType scan1( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

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
     *  @param[out]    perm       if not NULL, contains the permutation vector
     *  @param[out]    outValues  if not NULL, contains the sorted values
     *  @param[out]    inValues   array with the values to be sorted
     *  @param[in]     ascending  sort ascending (true) or descending (false)
     *  @param[in]     prefLoc    is the preferred context where computation should be done
     *
     *  Note: outValues = inValues[ perm ]
     *        i.e.: perm[i] contains the original position
     */
    template<typename ValueType>
    static void sort(
        hmemo::HArray<IndexType>* perm,
        hmemo::HArray<ValueType>* outValues,
        const hmemo::HArray<ValueType>& inValues,
        const bool ascending,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Sort an array of keys( index type ) and values in-place
     *
     *  @param[in,out] indexes    specify the keys used for sorting
     *  @param[in,out] values     sorted values
     *  @param[in]     ascending  sort ascending (true) or descending (false)
     *  @param[in]     prefLoc    is the preferred context where computation should be done
     *
     */
    template<typename ValueType>
    static void sortSparseEntries(
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<ValueType>& values,
        const bool ascending,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

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

    /** Determine the inverse of a permutation
     *  
     *  @param[in] perm is a permutation of the index values 0, .., n-1 
     *  @param[out] invPerm contains the inverse permutation 
     *  @param[in] prefLoc is the preferred context where computation should be done
     *  
     *  This routine throws an exception, if perm is not a valid permutation
     *
     *  \code
     *     perm    = { 0, 1, 3, 2, 5, 6, 7, 4 } ;
     *     invPerm = { 0, 1, 3, 2, 7, 4, 5, 6 } ;
     *     perm[invPerm] = { 0, 1, 2, 3, 4, 5, 6, 7 };
     *     invPerm[perm] = { 0, 1, 2, 3, 4, 5, 6, 7 };
     *  \endcode
     */
    static void inversePerm( 
        hmemo::HArray<IndexType>& invPerm,
        const hmemo::HArray<IndexType>& perm,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Sorting of an array where its subarrays are already sorted
     *
     *  @param[in,out] values array to be sorted
     *  @param[in,out] perm array where values reorded in the same way as values
     *  @param[in]     offsets array with offsets that specify the sorted subarrays
     *  @param[in]     ascending true for ascending sort, false for descending
     *  @param[in]     prefLoc context where merging should take place
     *
     *  Note: offset[0] == 0, offset[ offset.size() ] == array.size() must be valid
     *
     *  Each array[offset[i] : offset[i+1]] must already be sorted.
     */
    template<typename ValueType>
    static void mergeSort(
        hmemo::HArray<ValueType>& values,
        hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<IndexType>& offsets,
        bool ascending,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static void mergeSort(
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<IndexType>& offsets,
        bool ascending,
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
    static void setSequence(
        hmemo::HArray<ValueType>& array,
        ValueType startValue,
        ValueType inc,
        IndexType n,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Set an array with random values, untyped version.
     *
     *  @param[out] array    arbitray array, will contain random values of its type
     *  @param[in]  n        number of values, becomes size of array
     *  @param[in]  fillRate ratio of non-zero values
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */

    static void setRandom( hmemo::_HArray& array,
                           IndexType n,
                           float fillRate = 1.0f,
                           hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Set an array with random values, typed version
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

    static void buildSparseArray(
        hmemo::_HArray& sparseArray,
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::_HArray& denseArray,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build sparse indexes only, useful if sparseArray is not really needed */

    template<typename TargetType, typename SourceType>
    static void buildSparseArrayImpl(
        hmemo::HArray<TargetType>& sparseArray,
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::HArray<SourceType>& denseArray,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build sparse indexes only, useful if sparseArray is not really needed */

    template<typename ValueType>
    static void buildSparseIndexes(
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::HArray<ValueType>& denseArray,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Add new value in a sorted array, larger values are shifted up
     * 
     *  @param[in,out] array contains sorted values (ascending)
     *  @param[in]     value is the new entry to be inserted                          
     *  @param[in]     prefLoc  optional the context where the sparseIndexes are updated
     *  @returns       position in array where the entry has been added
     * 
     *  Note: the size of the array is increased by 1
     */
    template<typename ValueType>
    static IndexType insertSorted( 
        hmemo::HArray<ValueType>& array,
        const ValueType value,
        hmemo::ContextPtr = hmemo::ContextPtr() );

    /** Insert new value in an array at a certain pos, all other elements are shifted up
     * 
     *  @param[in,out] values
     *  @param[in]     pos is the position where the value is added
     *  @param[in]     val is the inserted value
     *  @param[in]     prefLoc  optional the context where the values are updated
     */
    template<typename ValueType>
    static void insertAtPos(
        hmemo::HArray<ValueType>& array,
        const IndexType pos,
        const ValueType value,
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

    static void buildDenseArray(
        hmemo::_HArray& denseArray,
        const IndexType denseN,
        const hmemo::_HArray& sparseArray,
        const hmemo::HArray<IndexType>& sparseIndexes,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Find an index in an array of sorted indexes 
     *
     *  @param[in] indexes is an array of sorted indexes
     *  @param[in] index is the value to be found in indexes array
     *  @return position of index in indexes if found, nIndex otherwise
     *
     *  \code
     *     HArray<IndexType> indexes ( 5, { 0, 5, 11, 18, 19 } );
     *     findPosInSortedIndexes( indexes, 5  ) -> 1
     *     findPosInSortedIndexes( indexes, 31 ) -> nIndex
     *  \endcode
     */
    static IndexType findPosInSortedIndexes( const hmemo::HArray<IndexType>& indexes, const IndexType pos );

    /** Add two sparse arrays 
     *
     *  @param[out] resultIndexes, resultValues for sparse result array
     *  @param[in]  indexes1, values1, zero1 stand for first sparse array
     *  @param[in]  alpha scaling factor for values of first array
     *  @param[in]  indexes2, values2, zero2 stand for second sparse array
     *  @param[in]  beta scaling factor for values of second array
     *  @param[in]  prefLoc is the context where operation should be done
     *
     *  Alias of any input array with the output array is not allowed, e.g. 
     *  resultValues and values1 must not be the same array.
     */
    template<typename ValueType>
    static void addSparse(
        hmemo::HArray<IndexType>& resultIndexes,
        hmemo::HArray<ValueType>& resultValues,
        const hmemo::HArray<IndexType>& indexes1,
        const hmemo::HArray<ValueType>& values1,
        const ValueType zero1,
        const ValueType alpha,
        const hmemo::HArray<IndexType>& indexes2,
        const hmemo::HArray<ValueType>& values2,
        const ValueType zero2,
        const ValueType beta,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise binary operation on sparse array: result[i] = op( x[i], y[i] )
     *
     *  @param[out] resultIndexes, resultValues for sparse result array
     *  @param[in]  indexes1, values1, zero1 stand for first sparse array
     *  @param[in]  indexes2, values2, zero2 stand for second sparse array
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     */
    template<typename ValueType>
    static void binaryOpSparse(
        hmemo::HArray<IndexType>& resultIndexes,
        hmemo::HArray<ValueType>& resultValues,
        const hmemo::HArray<IndexType>& indexes1,
        const hmemo::HArray<ValueType>& values1,
        const ValueType zero1,
        const hmemo::HArray<IndexType>& indexes2,
        const hmemo::HArray<ValueType>& values2,
        const ValueType zero2,
        const common::binary::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

private:

    template<typename ValueType>
    static void mergeSortOptional(
        hmemo::HArray<ValueType>& values,
        hmemo::HArray<IndexType>* perm,
        const hmemo::HArray<IndexType>& offsets,
        bool ascending,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    HArrayUtils();  // static class, no objects outside
    HArrayUtils( const HArrayUtils& );

};

} /* end namespace utilskernel */

} /* end namespace scai */

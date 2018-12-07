/**
 * @file HArrayUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
#include <scai/common/CompareOp.hpp>
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
    static void _assign(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename TargetValueType, typename SourceValueType>
    static void assign(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** 
     *  Assign to a resized array, so either target array is filled up or truncated.
     */
    template<typename TargetValueType, typename SourceValueType>
    static void assignResized(
        hmemo::HArray<TargetValueType>& target,
        const IndexType newSize,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read) of values with heterogeneous arrays.
     *
     *  target[i] op= source[indexes[i]]
     */
    static void _gather(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering of values in sparse array
     *
     *  This operation is the same as gather but here the source array is reprensented
     *  as a sparse version.
     *
     *  @param target will contain results, can be any type
     *  @param sourceZeroValue default value of source array if index does not appear in sourceNonZeroIndexes
     *  @param sourceNonZeroIndexes indexes of source array with non-zero values
     *  @param sourceNonZeroValues same size as sourceNonZeroIndexes, contains values of source array
     *  @param indexes same size as target, indexes of source array to gather
     *  @param op specifies binary op how to combine with existing values, can be COPY
     *  @param prefLoc preferred location where to execute the operation
     *
     *  \code
     *     for ( i = 0; i < indexes.size(); ++i )
     *       if ( there is j with sourceNonZeroIndexes[j] == indexes[i] )
     *         target[i] = target[i] op sourceNonZeroValues[j] ;
     *       else
     *         target[i] = target[i] op sourceZeroValue;
     *  \endcode
     *
     *  As sourceNonZeroIndexes is sorted (ascending), binary search can be applied to find j for i
     */
    template<typename SourceValueType>
    static void _sparseGather(
        hmemo::_HArray& target,
        const SourceValueType sourceZeroValue,
        const hmemo::HArray<SourceValueType>& sourceNonZeroValues,
        const hmemo::HArray<IndexType>& sourceNonZeroIndexes,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read of values) with HArrays, template typed version
     */
    template<typename TargetValueType, typename SourceValueType>
    static void gather(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Gathering (unstructured read of values) for sparse data, typed version
     */
    template<typename TargetValueType, typename SourceValueType>
    static void sparseGather(
        hmemo::HArray<TargetValueType>& target,
        const SourceValueType sourceZero,
        const hmemo::HArray<SourceValueType>& sourceNonZeroValues,
        const hmemo::HArray<IndexType>& sourceNonZeroIndexes,
        const hmemo::HArray<IndexType>& indexes,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /**
     *  @brief Scattering (unstructured write) of values with heterogeneous arrays.
     *
     *  target[index[i]] = source[i]
     */
    static void _scatter(
        hmemo::_HArray& target,
        const hmemo::HArray<IndexType>& index,
        const bool unique,
        const hmemo::_HArray& source,
        const common::BinaryOp op,
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
        const bool unique,
        const hmemo::HArray<SourceValueType>& source,
        const common::BinaryOp op,
        const hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** This method sets a single value in a heterogeneous array.
     *
     *  @param[in,out] target  Harray where a value to set
     *  @param[in]     index  position to set ( 0 <= index < target.size() )
     *  @param[in]     val    value to set
     *  @param[in]     op     specifies how to combine with current value
     *
     *  The value will be set at a valid context.
     */
    template<typename ValueType>
    static void setVal(
        hmemo::HArray<ValueType>& target,
        const IndexType index,
        const ValueType val,
        const common::BinaryOp op = common::BinaryOp::COPY );

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
     *
     *  Note: size of x and y must be equal if alpha != 0 and beta != 0
     *  Note: size of result is x.size() if alpha != 0, otherwise y.size()
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
     *  Any alias of the arrays x, y, result is supported. If no alias exists
     *  result will be allocated with the same size as x and y.
     */

    template<typename ValueType>
    static void arrayTimesArray(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** @brief Append an array to a given array 
     *
     *  @param[in,out] array1  is the array that will be extended
     *  @param[in]     array2  is the array that willl be appended
     *  @param[in]     context array1 should be valid here afterwards 
     */
    template<typename ValueType>
    static void appendArray(
        hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr context = hmemo::ContextPtr() );
    /*
     * Implementation of functions
     */
    static void _setArray(
        hmemo::_HArray& target,
        const hmemo::_HArray& source,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    template<typename TargetValueType, typename SourceValueType>
    static void setArray(
        hmemo::HArray<TargetValueType>& target,
        const hmemo::HArray<SourceValueType>& source,
        const common::BinaryOp op = common::BinaryOp::COPY,
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

    static void _setArraySection(
        hmemo::_HArray& target,
        const IndexType targetOffset,
        const IndexType targetStride,
        const hmemo::_HArray& source,
        const IndexType sourceOffset,
        const IndexType sourceStride,
        const IndexType n,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    /** Typed version for setting sectioned arrays */

    template<typename TargetValueType, typename SourceValueType>
    static void setArraySection(
        hmemo::HArray<TargetValueType>& target,
        const IndexType targetOffset,
        const IndexType targetStride,
        const hmemo::HArray<SourceValueType>& source,
        const IndexType sourceOffset,
        const IndexType sourceStride,
        const IndexType n,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    /** Fill a regular section of an array with a certain value */

    template<typename ValueType>
    static void fillArraySection(
        hmemo::HArray<ValueType>& value,
        const IndexType offset,
        const IndexType stride,
        const ValueType val,
        const IndexType n,
        const common::BinaryOp op = common::BinaryOp::COPY,
        hmemo::ContextPtr context = hmemo::ContextPtr() );

    template<typename ValueType>
    static void setScalar(
        hmemo::HArray<ValueType>& target,
        const ValueType value,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Computing the transpose of two-dimensional array.
     *
     *  @param[out] target result array, size will be n1 * n2
     *  @param[in]  n1 number of rows of 2D array
     *  @param[in]  n2 number of cols of 2D array
     *  @param[in]  source array with size n2 * n1
     *  @param[in]  conj   if true, compuate conjugate-transpose
     *  @param[in]  prefLoc location where result is computed.
     *
     *  Alias of target and source is suppported. On the host, 
     *  this operation is done in-place, on a device a temporary array
     *  is used.
     */
    template<typename ValueType>
    static void transpose(
        hmemo::HArray<ValueType>& target,
        const IndexType n1, 
        const IndexType n2, 
        const hmemo::HArray<ValueType>& source,
        const bool conj,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );
  
    /** Initialize an array with a certain size and a given value.
     *
     *  @param[out] array    is the array that will be allocated and set
     *  @param[in]  n        is the number of values
     *  @param[in]  val      is the value for each entry of the array
     *  @param[in]  prefLoc  is the location where the array is initialized, defaults to Host
     */
    template<typename ValueType>
    static void setSameValue(
        hmemo::HArray<ValueType>& array,
        const IndexType n,
        const ValueType val,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Provide a version of setSameValue that deals with untyped HArray. */

    template<typename ValueType>
    static void _setSameValue(
        hmemo::_HArray& target,
        const IndexType n,
        const ValueType val,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType getVal(
        const hmemo::HArray<ValueType>& array,
        const IndexType index );

    template<typename ValueType>
    static ValueType reduce(
        const hmemo::HArray<ValueType>& array,
        const common::BinaryOp redOp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static inline ValueType min(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static inline ValueType max(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static inline RealType<ValueType> maxNorm(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static inline ValueType sum(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType reduce2(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        const common::BinaryOp binOp,
        const common::BinaryOp redOp,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Functions that returns true if the element-wise comparison 
     *  of array elements returns true for all entries. 
     *
     *  @param[in] array1, array2 input arrays must have same size
     *  @param[in] compareOp specifies comparison operator to use
     *  @param[in] prefLoc optional the context where operation should take place
     */
    template<typename ValueType>
    static bool all(
        const hmemo::HArray<ValueType>& array1,
        const common::CompareOp compareOp,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Functions that returns true if the comparison of each array
     *  element with a scalar returns true.
     *
     *  @param[in] array input array
     *  @param[in] value is the scalar entry against which elements of array are compared
     *  @param[in] compareOp specifies comparison operator to use
     *  @param[in] prefLoc optional the context where operation should take place
     */
    template<typename ValueType>
    static bool allScalar(
        const hmemo::HArray<ValueType>& array,
        const common::CompareOp compareOp,
        const ValueType value,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static RealType<ValueType> l1Norm(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static RealType<ValueType> l2Norm(
        const hmemo::HArray<ValueType>& array,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static RealType<ValueType> maxDiffNorm(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static ValueType dotProduct(
        const hmemo::HArray<ValueType>& array1,
        const hmemo::HArray<ValueType>& array2,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise UnaryOp operation on array: result[i] = op( x[i] )
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
        const common::UnaryOp op,
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
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Elementwise binary operation on array: result[i] = op( x[i], y ), second arg is scalar
     *
     *  @param[out] result  output array
     *  @param[in]  x       input array
     *  @param[in]  y       input value
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  swapScalar if true scalar is first arg               
     *  @param[in]  prefLoc location where operation should be done if possible
     *
     *  Note: swapScalar does only matter if op is not commutative
     */
    template<typename ValueType>
    static void binaryOpScalar(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const ValueType y,
        const common::BinaryOp op,
        const bool swapScalar,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** More convenient interface for binaryOpScalar( swap = true ) */

    template<typename ValueType>
    static void compute(
        hmemo::HArray<ValueType>& result,
        const hmemo::HArray<ValueType>& x,
        const common::BinaryOp op,
        const ValueType y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() )
    {
        binaryOpScalar( result, x, y, op, false, prefLoc );
    }

    /** More convenient interface for binaryOpScalar( swap = false ) */

    template<typename ValueType>
    static void compute(
        hmemo::HArray<ValueType>& result,
        const ValueType x,
        const common::BinaryOp op,
        const hmemo::HArray<ValueType>& y,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() )
    {
        // y and x must be swapped
        binaryOpScalar( result, y, x, op, true, prefLoc );
    }

    /** Check for an index array whether all values are smaller than n */

    static bool validIndexes( 
        const hmemo::HArray<IndexType>& array, 
        const IndexType size, 
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Check whether values in array are sorted  
     *
     *  @param[in] array
     *  @param[in] op
     *  @param[in] prefLoc 
     *  @returns   true if a[i] op a[i+1] is true for all neighbored pairs
     */
    template<typename ValueType>
    static bool isSorted( 
        const hmemo::HArray<ValueType>& array, 
        const common::CompareOp op,
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
    static void bucketSortOffsets(
        hmemo::HArray<IndexType>& offsets,
        hmemo::HArray<IndexType>& perm,
        const hmemo::HArray<BucketType>& array,
        const BucketType nb,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename BucketType>
    static void bucketSortSizes(
        hmemo::HArray<IndexType>& sizes,
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

    /** Eliminate duplicate elements in a sorted array 
     *
     *  @param[in,out] indexes is the sorted array of indexes that might contain double values
     *  @param[in,out] values are values that will be combined
     *  @param[in] op  specifies how to combine values with same index pos
     *
     *  /code
     *    in:    indexes[] = { 0, 1, 5, 7, 7, 9, 9 };    
     *    in:    values [] = { 0, 1, 2, 3, 4, 5, 6 };
     *    elimDoubles ( indexes, values, BinaryOp::COPY )
     *    out:   indexes[] = { 0, 1, 5, 7, 9 };    
     *    out:   values [] = { 0, 1, 2, 4, 6 };
     *    elimDoubles ( indexes, values, BinaryOp::ADD )
     *    out:   indexes[] = { 0, 1, 5, 7, 9 };    
     *    out:   values [] = { 0, 1, 2, 7, 11 };
     *  /endcode
     */
    template<typename ValueType>
    static void elimDoubles(
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<ValueType>& values,
        const common::BinaryOp op );

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
     *  @param[out] array    arbitray array, is filled with random values of its type
     *  @param[in]  bound    random values between 0 and bound
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */
    static void setRandom( hmemo::_HArray& array,
                           IndexType bound,
                           hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Replace entries of an array with random values with a certain probability for each element
     *
     *  @param[out] array    arbitray array, is filled with random values of its type
     *  @param[in]  fillRate probability whether one array element is filled or remains unchanged
     *  @param[in]  bound    random values between 0 and bound
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */
    static void setSparseRandom( hmemo::_HArray& array,
                                 float fillRate,
                                 IndexType bound,
                                 hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Set an array with random values, typed version
     *
     *  @param[out] array    will contain random values of its type
     *  @param[in]  bound    random values between 0 and bound
     *  @param[in]  fillRate probability whether one array element is filled or remains unchanged
     *  @param[in]  prefLoc  optional the context where random numbers should be drawn
     */
    template<typename ValueType>
    static void fillRandom( hmemo::HArray<ValueType>& array,
                            IndexType bound,
                            float fillRate,
                            hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Create an array of sparse indexes where an index value appears only with a certain probability.
     *
     *  @param[out] array       with sparse indexes, maximal size is n
     *  @param[in]  n           is upper bound of range, i.e. indexes are 0 to n-1
     *  @param[in]  probability whether an index appears or not, 0 for never, 1 for always.
     */
    static void randomSparseIndexes( hmemo::HArray<IndexType>& array, const IndexType n, const float probability );

    /** Build sparse array from dense array, needed for conversion DenseVector -> SparseVector 
     *
     *  @param[out] sparseArray contains the non-zero values
     *  @param[out] sparseIndexes contains the indexes of sparseArray in original array
     *  @param[in]  denseArray    contains the array with all values
     *  @param[in]  zeroValue     is the value that is considered as zero in denseArray
     *  @param[in]  prefLoc       is the preferred location where operation is done
     * 
     *  \code
     *    buildSparseArray( sparseArray, sparseIndexes, HArray<double>( { 0, 2, 1, 2, 1 } ), 1 )
     *    -> sparseArray = HArray<double> ( { 0, 2, 2 } ); sparseIndexes = HArray<indexes> ( { 0, 1, 3 } )
     *    buildSparseArray( sparseArray, sparseIndexes, HArray<double>( { 0, 2, 1, 2, 1 } ), 2 )
     *    -> sparseArray = HArray<double> ( { 0, 1, 1 } ); sparseIndexes = HArray<indexes> ( { 0, 2, 4 } )
     *  \endcode
     */
    template<typename TargetType, typename SourceType>
    static void buildSparseArray(
        hmemo::HArray<TargetType>& sparseArray,
        hmemo::HArray<IndexType>& sparseIndexes,
        const hmemo::HArray<SourceType>& denseArray,
        const SourceType zeroValue,
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
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Insert new value in an array at a certain pos, all other elements are shifted up
     * 
     *  @param[in,out] array is the array where to insert the element
     *  @param[in]     pos is the position where the value is added
     *  @param[in]     value is the inserted value
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
     *  @param[out]  denseArray    will contain the dense data (input its only its size)
     *  @param[in]   denseN        is the size of the dense array
     *  @param[in]   sparseArray   contains non-zero values
     *  @param[in]   sparseIndexes are the positions of the non-zero values
     *  @param[in]   zero          is the default value for all positions not specifiedin sparseIndexes
     *  @param[in]   prefLoc       is the context where operation should be done
     *
     *  Note: sparseIndexes must contain only indexes between 0 and denseN - 1
     */
    template<typename ValueType>
    static void buildDenseArray(
        hmemo::HArray<ValueType>& denseArray,
        const IndexType denseN,
        const hmemo::HArray<ValueType>& sparseArray,
        const hmemo::HArray<IndexType>& sparseIndexes,
        const ValueType zero,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Find an index in an array of sorted indexes 
     *
     *  @param[in] indexes is an array of sorted indexes
     *  @param[in] index is the value to be found in indexes array
     *  @return position of index in indexes if found, invalidIndex otherwise
     *
     *  \code
     *     HArray<IndexType> indexes ( 5, { 0, 5, 11, 18, 19 } );
     *     findPosInSortedIndexes( indexes, 5  ) -> 1
     *     findPosInSortedIndexes( indexes, 31 ) -> invalidIndex
     *  \endcode
     */
    static IndexType findPosInSortedIndexes( const hmemo::HArray<IndexType>& indexes, const IndexType index );

    /** Find many indexes in an array of sorted indexes */

    static void findPosInSortedIndexesV( hmemo::HArray<IndexType>& outPos,
                                         const hmemo::HArray<IndexType>& indexes,
                                         const hmemo::HArray<IndexType> inPos );

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
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    template<typename ValueType>
    static IndexType allSparse(
        bool& allFlag,
        const hmemo::HArray<IndexType>& indexes1,
        const hmemo::HArray<ValueType>& values1,
        const ValueType zero1,
        const hmemo::HArray<IndexType>& indexes2,
        const hmemo::HArray<ValueType>& values2,
        const ValueType zero2,
        const common::CompareOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Merge non-zero entries of sparse array
     *
     *  @param[out] resultIndexes, resultValues for sparse result array
     *  @param[in]  indexes1, values1 stand for first sparse array
     *  @param[in]  indexes2, values2  stand for second sparse array
     *  @param[in]  op      specifies operation to apply on input values
     *  @param[in]  prefLoc location where operation should be done if possible
     *
     *  Example can be found in HArrayUtilsTest
     */
    template<typename ValueType>
    static void mergeSparse(
        hmemo::HArray<IndexType>& resultIndexes,
        hmemo::HArray<ValueType>& resultValues,
        const hmemo::HArray<IndexType>& indexes1,
        const hmemo::HArray<ValueType>& values1,
        const hmemo::HArray<IndexType>& indexes2,
        const hmemo::HArray<ValueType>& values2,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Select real or imaginary part for an array of (complex) values 
     *
     *  @param[out] realValues gets same size as complexValues, will contain real or imaginary parts
     *  @param[in]  complexValues input array with complex values
     *  @param[in]  kind specifies whether to select real or imaginary part 
     *  @param[in]  prefLoc location where selection should be done if possible
     */
    template<typename ValueType>
    static void selectComplexPart(
        hmemo::HArray<RealType<ValueType> >& realValues,
        const hmemo::HArray<ValueType>& complexValues,
        const common::ComplexPart kind,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Build complex array by combining two real array, one for real and one for imag part
     *
     *  @param[out] complexValues gets same size as realValues and imagValues
     *  @param[in]  realValues    input array with real parts for complex values
     *  @param[in]  imagValues    input array with imag parts for complex values
     *  @param[in]  prefLoc location where selection should be done if possible
     */
    template<typename ValueType>
    static void buildComplex(
        hmemo::HArray<ValueType>& complexValues,
        const hmemo::HArray<RealType<ValueType> >& realValues,
        const hmemo::HArray<RealType<ValueType> >& imagValues,
        hmemo::ContextPtr prefLoc = hmemo::ContextPtr() );

    /** Merge two arrays and retain index maps from the original arrays into the result array.
     *
     * Merge the two *sorted* arrays `x` and `y` into the array `result`, overwriting any data in `result`,
     * and retain two separate maps from the indexes of `x` and `y` into the indexes of `result`.
     *
     * The input arrays x and y must be sorted according to the specified comparator, and the
     * result will also be sorted according to the same comparator.
     *
     *  @param[out] result     the resulting array after merging. Any data in the array will be overwritten. Size is x.size() + y.size().
     *  @param[out] xMap       an array of indices such that result[xMap[i]] = x[i] for i = 0, .., x.size() - 1.
     *  @param[out] yMap       an array of indices such that result[yMap[i]] = y[i] for i = 0, .., y.size() - 1.
     *  @param[in]  x          a *sorted* array.
     *  @param[in]  y          a *sorted* array.
     *  @param[in]  comparator the operator to use for the comparison, i.e. LE for ascending or GE for descending sort
     *  @param[in]  prefLoc    preferred context for computations.
     */
    template <typename ValueType>
    static void mergeAndMap(
        hmemo::HArray<ValueType>& result,
        hmemo::HArray<IndexType>& xMap,
        hmemo::HArray<IndexType>& yMap,
        const hmemo::HArray<ValueType>& x,
        const hmemo::HArray<ValueType>& y,
        const common::CompareOp comparator = common::CompareOp::LE,
        hmemo::ContextPtr prefLoc = hmemo::Context::getContextPtr() );

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

/* ============================================================================= */

template<typename ValueType>
ValueType HArrayUtils::min(
    const hmemo::HArray<ValueType>& values,
    hmemo::ContextPtr prefLoc )
{
    return HArrayUtils::reduce( values, common::BinaryOp::MIN, prefLoc );
}

template<typename ValueType>
ValueType HArrayUtils::max(
    const hmemo::HArray<ValueType>& values,
    hmemo::ContextPtr prefLoc )
{
    return HArrayUtils::reduce( values, common::BinaryOp::MAX, prefLoc );
}

template<typename ValueType>
ValueType HArrayUtils::sum(
    const hmemo::HArray<ValueType>& values,
    hmemo::ContextPtr prefLoc )
{
    return HArrayUtils::reduce( values, common::BinaryOp::ADD, prefLoc );
}

template<typename ValueType>
RealType<ValueType> HArrayUtils::maxNorm(
    const hmemo::HArray<ValueType>& values,
    hmemo::ContextPtr prefLoc )
{
    return HArrayUtils::reduce( values, common::BinaryOp::ABS_MAX, prefLoc );
}

} /* end namespace utilskernel */

} /* end namespace scai */

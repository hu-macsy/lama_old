/**
 * @file UtilKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA utilities on heterogeneous arrays provided as kernels.
 * @author Thomas Brandes
 * @date 03.04.2013
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/CompareOp.hpp>
#include <scai/common/UnaryOp.hpp>

namespace scai
{

/** Namespace for utilities on heterogeneous arrays (HArray) */

namespace utilskernel
{

/** Structure just to group traits for all Utils kernels.
 *
 *  This struct does not contain any data at all.
 *  Therefore it could also be a namespace but it is more convenient as
 *  each trait must always be used qualified: UtilKernelTrait::utiliy
 */

struct UtilKernelTrait
{
    /** @brief Trait for register kernel function validIndexes */

    struct validIndexes
    {
        /** Check that all values in array are in a certain range.
         *
         *  @param array is an array of index values
         *  @param n is the size of array
         *  @param size specifies the range in which array values must fit
         *  @return true if \f$ 0 \le array[i] < size \forall i = 0, ..., n-1\f$
         */

        typedef bool ( *FuncType )( const IndexType array[], const IndexType n, const IndexType size );
        static const char* getId()
        {
            return "Util.validIndexes";
        }
    };

    /** @brief Trait for register kernel function reduce that reduces elements of an array
     *
     *  @tparam ValueType specifies the value type used in the reduction.
     */
    template <typename ValueType>
    struct reduce
    {
        /** @brief reduce op for n contiguously stored valuels
         *
         *  @param[in] array is an array of values
         *  @param[in] n is the size of array
         *  @param[in] zero  is the zero element used for reduction
         *  @param[in] op is the binary reduction operator ( ADD for sum, MIN for minval, MAX for maxval, ...)
         *  @return reduced value corresponding to the reduction operator
         */

        typedef ValueType ( *FuncType ) ( const ValueType array[],
                                          const IndexType n,
                                          const ValueType zero,
                                          const common::BinaryOp op );
        static const char* getId()
        {
            return "Util.reduce";
        }
    };

    /** @brief Trait for register kernel function reduce that transforms/reduces elements of two array
     *
     *  @tparam ValueType specifies the value type used in the reduction.
     */
    template <typename ValueType>
    struct reduce2
    {
        /** @brief reduce combined values of two arays
         *
         *  @param[in] array1 is first array of values
         *  @param[in] array2 is second array of values
         *  @param[in] n is the size of arrays
         *  @param[in] binop is the binary operator applied elementwise on array1 and array2
         *  @param[in] zero  is the zero element used for reduction
         *  @param[in] redop is the binary reduction operator ( ADD for sum, MIN for minval, MAX for maxval, ...)
         *  @return reduced value corresponding to the reduction operator
         */

        typedef ValueType ( *FuncType ) ( const ValueType array1[],
                                          const ValueType array2[],
                                          const IndexType n,
                                          const common::BinaryOp binop,
                                          const ValueType zero,
                                          const common::BinaryOp redop );
        static const char* getId()
        {
            return "Util.reduce2";
        }
    };

    template <typename ValueType>
    struct allCompare
    {
        typedef bool ( *FuncType ) ( const ValueType array1[],
                                     const ValueType array2[],
                                     const IndexType n,
                                     const common::CompareOp op );
        static const char* getId()
        {
            return "Util.allCompare";
        }
    };

    template <typename ValueType>
    struct allCompareScalar
    {
        typedef bool ( *FuncType ) ( const ValueType array[],
                                     const ValueType scalar,
                                     const IndexType n,
                                     const common::CompareOp op );
        static const char* getId()
        {
            return "Util.allCompareScalar";
        }
    };

    template <typename ValueType>
    struct isSorted
    {
        /** @brief Predicate that tests whether a sequene is sorted.
         *
         *  @param[in] array values to be checked
         *  @param[in] n number of values to check
         *  @param[in] op specifies comparison operator that must hold for each neighbored pair 
         *  @returns true iff \f$ a[i] \le a[i+1] \f$ (ascending=true) or \f$ a[i] \ge a[i+1] \f$ (ascending=false)
         */

        typedef bool ( *FuncType ) ( const ValueType array[], const IndexType n, const common::CompareOp op );
        static const char* getId()
        {
            return "Util.isSorted";
        }
    };

    /** @brief Structure with function pointer type defintions for setter methods.
     *
     *  @tparam ValueType specifies the value type used in the set operations.
     */

    template<typename ValueType>
    struct setVal
    {
        /** Set all elements of a contiguous array with a value.
         *  A binary operator like ADD, MULT can be used to combine the new value with the old value.
         */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n, const ValueType val, const common::BinaryOp op );
        static const char* getId()
        {
            return "Util.setVal";
        }
    };

    template<typename ValueType>
    struct scaleVectorAddScalar
    {
        /** Calculates array1 = alpha * array2 + beta (elementwise) */

        typedef void ( *FuncType ) ( ValueType array1[], const ValueType array2[], const IndexType n, const ValueType alpha, const ValueType beta );
        static const char* getId()
        {
            return "Util.scaleVectorAddScalar";
        }
    };

    template<typename ValueType>
    struct setOrder
    {
        /** Set all elements of a contiguous array with its order number 0, 1, 2, ... */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n );
        static const char* getId()
        {
            return "Util.setOrder";
        }
    };

    template<typename ValueType>
    struct setSequence
    {
        /** Set all elements of a contiguous array with its order number 0, 1, 2, ... */

        typedef void ( *FuncType ) ( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );
        static const char* getId()
        {
            return "Util.setSequence";
        }
    };

    template<typename ValueType>
    struct getValue
    {
        typedef ValueType ( *FuncType ) ( const ValueType* array, const IndexType i );
        static const char* getId()
        {
            return "Util.getValue";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct set
    {
        /** Set out[i] _op= in[i],  0 <= i < n , op = +, -, *, /, min, max, ... */

        typedef void ( *FuncType ) ( TargetValueType out[], const SourceValueType in[], const IndexType n, const common::BinaryOp op );
        static const char* getId()
        {
            return "Util.set";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct setSection
    {
        /** Set out[i * inc_out] _op= in[i * inc_in],  0 <= i < n , op = +, -, *, /, min, max, ... */

        typedef void ( *FuncType ) ( TargetValueType out[], const IndexType inc_out,
                                     const SourceValueType in[], const IndexType inc_in,
                                     const IndexType n, const common::BinaryOp op );
        static const char* getId()
        {
            return "Util.setSection";
        }
    };

    template<typename ValueType>
    struct fillSection
    {
        /** Set out[i * inc] _op= val,  0 <= i < n , op = +, -, *, /, min, max, ... */

        typedef void ( *FuncType ) ( ValueType out[], const IndexType inc,
                                     const ValueType val,
                                     const IndexType n, const common::BinaryOp op );
        static const char* getId()
        {
            return "Util.fillSection";
        }
    };

    template<typename ValueType>
    struct unaryOp
    {
        /** Apply UnaryOp op sin/cos/sqrt/... function elementwise on vector
         *  This routine can also be used for aliased arrays, i.e. in == out
         *  This method can only be used for numeric types, not for IndexType
         */
        typedef void ( *FuncType ) (
            ValueType out[],
            const ValueType in[],
            const IndexType n,
            const common::UnaryOp op );

        static const char* getId()
        {
            return "Util.unaryOp";
        }
    };

    template<typename ValueType>
    struct binaryOp
    {
        /** Apply binary op ADD, MULT, POW, COPY_SIGN, ... on array with one given numeric type.
         *  This routine can also be used for aliased arrays, i.e. in1 == out or in2 == out
         *  This method can only be used for numeric types, not for IndexType
         *
         *  @param[out] out array with output values
         *  @param[in]  in1 array with input values
         *  @param[in]  in2 array with input values
         *  @param[in]  n   number of elements
         *  @param[in]  op  binary operation to be applied
         */
        typedef void ( *FuncType ) (
            ValueType out[],
            const ValueType in1[],
            const ValueType in2[],
            const IndexType n,
            const common::BinaryOp op );

        static const char* getId()
        {
            return "Util.binaryOp";
        }
    };

    template<typename ValueType>
    struct binaryOpScalar
    {
        /** Same as binaryOp but one operand is only a scalar
         *
         *  This operation is only available for numeric types, not for IndexType
         *
         *  @param[out] out        array with output values
         *  @param[in]  in         array with values for second argument
         *  @param[in]  value      scalar value as first argument
         *  @param[in]  n          number of elements
         *  @param[in]  op         binary operation to be applied
         *  @param[in]  swapScalar if true value becomes first argument of binary op
         */
        typedef void ( *FuncType ) (
            ValueType out[],
            const ValueType in[],
            const ValueType value,
            const IndexType n,
            const common::BinaryOp op,
            const bool swapScalar );

        static const char* getId()
        {
            return "Util.binaryOpScalar";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct setGather
    {
        /** Set out[i] op = in[ indexes[i] ],  \f$0 \le i < n\f$ */

        typedef void ( *FuncType ) (
            TargetValueType out[],
            const SourceValueType in[],
            const IndexType indexes[],
            const common::BinaryOp op,
            const IndexType n );

        static const char* getId()
        {
            return "Util.setGather";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct setGatherSparse
    {
        /** setGather where in is a sparse array */

        typedef void ( *FuncType ) (
            TargetValueType out[],
            const SourceValueType inZeroValue,
            const SourceValueType inNonZeroValues[],
            const IndexType inNonZeroIndexes[],
            const IndexType nnz,
            const IndexType indexes[],
            const common::BinaryOp op,
            const IndexType n );

        static const char* getId()
        {
            return "Util.setGatherSparse";
        }
    };

    template<typename ValueType>
    struct scatterVal
    {
        /** @brief Setting one value at multiple positions in an array.
         *
         *  @param[in,out] out is the array in which values will be written
         *  @param[in]     indexes are the positions the value is written
         *  @param[in]     value is the value written in the out array
         *  @param[in]     n is the number of values
         *
         *  Note: Not all values might be set in 'out'.
         *
         *  out[ indexes[i] ] = value, i = 0, ..., n-1
         */

        typedef void ( *FuncType ) (
            ValueType out[],
            const IndexType indexes[],
            const ValueType value,
            const IndexType n );

        static const char* getId()
        {
            return "Util.scatterVal";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct setScatter
    {
        /** @brief Indirect set of arrays also known as scatter.
         *
         *  @param[in,out] out is the array in which values will be inserted
         *  @param[in]     indexes are the positions where values are written
         *  @param[in]     unique if true no index appears twice
         *  @param[in]     in is the array with the output values.
         *  @param[in]     op specifies how the set element is combined with available element
         *  @param[in]     n is the number of values
         *
         *  Note: If indexes are unique, the operation is done data parallel;
         *        if they are not unique, a parallel execution is only possible if
         *        the binary operation can be done via atomic updates.
         *
         *  out[ indexes[i] ] = in [i] , i = 0, ..., n-1   for op == recution::COPY
         *  out[ indexes[i] ] += in [i] , i = 0, ..., n-1   for op == recution::ADD
         *  out[ indexes[i] ] *= in [i] , i = 0, ..., n-1   for op == recution::MULT
         */

        typedef void ( *FuncType ) (
            TargetValueType out[],
            const IndexType indexes[],
            const bool unique,
            const SourceValueType in[],
            const common::BinaryOp op,
            const IndexType n );

        static const char* getId()
        {
            return "Util.setScatter";
        }
    };

    template<typename ValueType>
    struct scan
    {
        /** This method computes runnings sums of values
         *
         *  @param[in,out] array contains  values and later the running sums
         *  @param[in]    n is the number of values 
         *  @param[in]    first is added to all elements
         *  @param[in]    exclusive if false a[i] = a[i] op a[i-1], else a[i] = a[i-1] op a[i-2] op ... op first
         *  @param[in]    append if true the result value is added at the end of the array
         *  @returns      the total sum of values
         *
         *  \code
         *    array  :    3   7   8   4   2   
         *    array  :    3  10  18  22  24    -> returns 24, exclusive = false
         *    array  :    0   3  10  18  22    -> returns 24, exclusive = true
         *  \endcode
         *
         *  \code
         *     array[i] = array'[i-1] + array'[i-2] + ... + array'[0] + first   exclusive = true
         *     array[i] = array'[i] + array'[i-1] + ... + array'[0] + first     exclusive = false
         *  \endcode
         *
         *  Important: array must have n + 1 allocated entries if append is true
         *
         */

        typedef ValueType ( *FuncType ) ( ValueType array[], const IndexType n, const ValueType first, bool exclusive, bool append );

        static const char* getId()
        {
            return "Util.scan";
        }
    };

    template<typename ValueType>
    struct unscan
    {
        /** This method computes differences of values, array[i] = array[i+1] - array[i]
         *
         *  @param[in,out] array contains  values and later the differences
         *  @param[in]     n is the number of values, array must contain at least this number of vales
         *  @returns       0
         *
         *  \code
         *    array  :    0   3  10  18  22  24
         *    array  :    3   7   8   4   2  x   -> returns 0
         *  \endcode
         *
         *  Important: sizes must have numRows + 1 allocated entries.
         *
         */

        typedef ValueType ( *FuncType ) ( ValueType array[], const IndexType n );

        static const char* getId()
        {
            return "Util.unscan";
        }
    };

    template<typename ValueType>
    struct sort
    {
        /** Stable sorting of values
         *
         *  @param[out]  perm      contains the positions of the sorted values in the input array inValues, optional
         *  @param[out]  outValues array with sorted values, optional
         *  @param[in]   inValues  array with values to be sorted
         *  @param[in]   n is the number of values to be sorted
         *  @param[in]   ascending if true sort in ascending order, descending otherwise
         *
         *  \code
         *         inValues =   1  4   1  8  5  7
         *
         *        outValues =   8   7   5   4   1   1
         *           perm   =   3   5   4   1   0   2
         *
         *         outValues = inValues[ perm ]
         *  \endcode
         *
         *  - perm == NULL specifies that this optional argument is not present
         *  - outValues == NULL specifies that this optional argument is not present
         *  - at least perm or outValues must be available
         *  - inValues == outValues is supported, this alias implies sorting in-place
         *
         *  Note: if the values to be sorted are int values in a certain range like 0, ..., nb - 1,
         *        sorting might be more efficient with countBuckets, sortInBuckets
         */
        typedef void ( *FuncType ) (
            IndexType perm[],
            ValueType outValues[],
            const ValueType inValues[],
            const IndexType n,
            const bool ascending );

        static const char* getId()
        {
            return "Utils.sort";
        }
    };

    template<typename ValueType>
    struct sortInPlace
    {
        /** Sorting of values
         *
         *  @param[in,out] indexes
         *  @param[in,out] values
         *  @param[in]   n is the number of values to be sorted
         *  @param[in]   ascending if true sort in ascending order, descending otherwise
         *
         *  \code
         *         indexes  =     8   3    1   5
         *         values   =   1.0  2.0  3.0  4.0
         *
         *         indexes  =     1   3    5    8
         *         values   =   3.0  2.0  4.0  1.0
         *  \endcode
         */
        typedef void ( *FuncType ) (
            IndexType indexes[],
            ValueType values[],
            const IndexType n,
            const bool ascending );

        static const char* getId()
        {
            return "Utils.sortInPlace";
        }
    };

    template<typename BucketType>
    struct countBuckets
    {
        /** Count bucket sizes for values mapped to buckets
         *
         *  @param[in] nBuckets  number of buckets
         *  @param[in] n number of values to sort in buckets
         *  @param[in] bucketMap array with n entries, bucketMap[i] is bucket for entry i
         *  @param[out] bucketSizes array with nBuckets entries, bucketSizes[j] is number of elements mapped to this bucket
         *
         *  \code
         *           bucketMap [n=10]           =  {  0  1  2 0  1  2  1  2  0  0 }
         *           bucketSizes [nBuckets = 3] =  {  4  3  3 }
         *  \endcode
         *
         *  Note: sum( bucketSizes ) = n implies that all entries in bucketMap were legal buckets between 0 and nBucket-1
         *        This routine does not throw an exception for illegal entries
         */
        typedef void ( *FuncType ) ( IndexType bucketSizes[], const BucketType nBuckets, const BucketType bucketMap[], const IndexType n );

        static const char* getId()
        {
            return "Utils.countBuckets";
        }
    };

    template<typename BucketType>
    struct sortInBuckets
    {
        /** Resort indexes 0, ..., n-1 according to their mapping to buckets
         *
         *  @param[in] nBuckets  number of buckets
         *  @param[in] n number of values to sort in buckets
         *  @param[in] bucketMap array with n entries, bucketMap[i] is bucket for entry i
         *  @param[in] offsets array with nBuckets + 1 entries, is running sum
         *
         *  \code
         *           bucketMap [n=10]               =  {  0  1  2  0  1  2  1  2  0  0 }
         *           bucketSizes [nBuckets = 3]     =  {  4  3  3 }
         *           offsets     [nBuckets + 1 = 4] =  {  0           4        7        10  }
         *           sortedIndexes [n=10]           =  {  0  3  8  9  1  4  6  2  5  7 }
         *  \endcode
         *
         */

        typedef void ( *FuncType )( IndexType sortedIndexes[],
                                    IndexType offsets[],
                                    const BucketType nBuckets,
                                    const BucketType bucketMap[],
                                    const IndexType n );

        static const char* getId()
        {
            return "Utils.sortInBuckets";
        }
    };

    struct setInversePerm
    {
        /** Compute the inverse permutation for a given permutation.
         *
         *  inversePerm [ perm [i] ] == i , 0 <= i < n
         *
         *  @param[out] inversePerm, size = n, will contain the inverse permutation
         *  @param[in] perm, size = n, is input permuation of 0, ..., n-1
         *  @param[in] n specifies the size of perm and inversePerm
         *
         *  /code
         *       perm      2  5  1  4  6  3  0
         *     inperm      6  2  0  5  3  1  4
         *  /endcode
         */

        typedef void ( *FuncType ) ( IndexType inversePerm[],
                                     const IndexType perm[],
                                     const IndexType n );

        static const char* getId()
        {
            return "Utils.setInversePerm";
        }
    };
};

} /* end namespace utilskernel */

} /* end namespace scai */

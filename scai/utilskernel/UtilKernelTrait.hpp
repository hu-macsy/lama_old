/**
 * @file UtilKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA utilities provided as kernels.
 * @author Thomas Brandes
 * @date 03.04.2013
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/ReductionOp.hpp>

namespace scai
{

/** Namespace for utilities on heterogeneous arrays (HArray) and derived class LArray */

namespace utilskernel
{

/** Structure with traits for all Utils kernels. */

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
         *  @param[in] op is the reduction operator ( ADD for sum, MIN for minval, MAX for maxval, ...)
         *  @return reduced value corresponding to the reduction operator
         */

        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n, const reduction::ReductionOp op );
        static const char* getId()
        {
            return "Util.reduce";
        }
    };

    template <typename ValueType>
    struct absMaxDiffVal
    {
        /** @brief Building absolute maximum of element-wise difference of vector elements.
         *
         *  @param array1i[in] first array
         *  @param array2i[in] second array
         *  @param n           size of array1 and array2
         *  @returns           max( abs( array1[i] - array2[i] ) ), \f$ 0 \le i < n \f$
         *
         *  Function is helpful to compute maximum norm for vectors and matrices
         */

        typedef ValueType ( *FuncType ) ( const ValueType array1[], const ValueType array2[], const IndexType n );
        static const char* getId()
        {
            return "Util.absMaxDiffVal";
        }
    };

    template <typename ValueType>
    struct copysign
    {
        /** @brief Composes a floating point value with the magnitude of x and the sign of y
         *
         *  @param x[out] result array
         *  @param x[in]  first array giving the magnitude of the result vector
         *  @param y[in]  second array giving the sign of the result vector
         *  @param n      size of the arrays
         *
         *  Function is helpful to compute maximum norm for vectors and matrices
         */

        typedef void ( *FuncType ) ( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n );
        static const char* getId()
        {
            return "Util.copysign";
        }
    };

    template <typename ValueType>
    struct isSorted
    {
        /** @brief Predicate that tests whether a sequene is sorted.
         *
         *  @param[in] array values to be checked
         *  @param[in] n number of values to check
         *  @param[in] ascending if true check for ascending order, otherwise for descending
         */

        typedef bool ( *FuncType ) ( const ValueType array[], const IndexType n, bool ascending );
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
         *  A reduction operator like ADD, MULT can be used to combine the new value with the old value.
         */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n, const ValueType val, const reduction::ReductionOp op );
        static const char* getId()
        {
            return "Util.setVal";
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

    template<typename ValueType1, typename ValueType2>
    struct set
    {
        /** Set out[i] _op= in[i],  0 <= i < n , op = +, -, *, /, min, max, ... */

        typedef void ( *FuncType ) ( ValueType1 out[], const ValueType2 in[], const IndexType n, const reduction::ReductionOp op );
        static const char* getId()
        {
            return "Util.set";
        }
    };

    template<typename ValueType1, typename ValueType2>
    struct setSection
    {
        /** Set out[i + inc_out] _op= in[i + inc_in],  0 <= i < n , op = +, -, *, /, min, max, ... */

        typedef void ( *FuncType ) ( ValueType1 out[], const IndexType inc_out,
                                     const ValueType2 in[], const IndexType inc_in,
                                     const IndexType n, const reduction::ReductionOp op );
        static const char* getId()
        {
            return "Util.setSection";
        }
    };

    template<typename ValueType>
    struct execElementwise
    {
        /** Execute sin/cos/sqrt/... function elementwise on vector */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n, const elementwise::ElementwiseOp op );
        static const char* getId()
        {
            return "Util.execElementwise";
        }
    };

    template <typename ValueType>
    struct vectorScale
    {
        /** @brief Building absolute maximum of element-wise difference of vector elements.
         *
         *  @param result[out] result array
         *  @param x[in]       first input array
         *  @param y[in]       second input array
         *  @param n           size of all three arrays
         *  @returns           result[i] = x[i] * y[i], \f$ 0 \le i < n \f$
         *
         *  Function is helpful to compute maximum norm for vectors and matrices
         */

        typedef void ( *FuncType ) ( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n );
        static const char* getId()
        {
            return "Util.vectorScale";
        }
    };

    template<typename ValueType1, typename ValueType2>
    struct setScale
    {
        /** @brief scaled array assignment, out = in * value
         *
         *  Set out[i] = scale * in[i],  0 <= i < n
         *
         *  @param[out]     outValues  is the output array
         *  @param[in]      scaleValue scaling factor
         *  @param[in]      inValues   is the array with entries to scale
         *  @param[in]      n          is the number of entries
         */
        typedef void ( *FuncType ) (
            ValueType1 outValues[],
            const ValueType1 scaleValue,
            const ValueType2 inValues[],
            const IndexType n );

        static const char* getId()
        {
            return "Util.setScale";
        }
    };

    template<typename ValueType1, typename ValueType2>
    struct setGather
    {
        /** Set out[i] = in[ indexes[i] ],  \f$0 \le i < n\f$ */

        typedef void ( *FuncType ) (
            ValueType1 out[],
            const ValueType2 in[],
            const IndexType indexes[],
            const IndexType n );

        static const char* getId()
        {
            return "Util.setGather";
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

    template<typename ValueType1, typename ValueType2>
    struct setScatter
    {
        /** @brief Indirect set of arrays also known as scatter.
         *
         *  @param[in,out] out is the array in which values will be inserted
         *  @param[in]     indexes are the positions where values are written
         *  @param[in]     in is the array with the output values.
         *  @param[in]     op specifies how the set element is combined with available element
         *  @param[in]     n is the number of values
         *
         *  Note: Not all values might be set in 'out'. There should be no double
         *        values in indexes as this might result in non-ambiguous results
         *        by a parallel execution.
         *
         *  out[ indexes[i] ] = in [i] , i = 0, ..., n-1   for op == recution::COPY
         *  out[ indexes[i] ] += in [i] , i = 0, ..., n-1   for op == recution::ADD
         *  out[ indexes[i] ] *= in [i] , i = 0, ..., n-1   for op == recution::MULT
         */

        typedef void ( *FuncType ) (
            ValueType1 out[],
            const IndexType indexes[],
            const ValueType2 in[],
            const reduction::ReductionOp op,
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
         *  @param[in]    n is the number of values, array must contain one additional value
         *  @returns      the total sum of values
         *
         *  \code
         *    array  :    3   7   8   4   2  x
         *    array  :    0   3  10  18  22  24      -> returns 24
         *  \endcode
         *
         *  Important: sizes must have numRows + 1 allocated entries.
         *
         */

        typedef ValueType ( *FuncType ) ( ValueType array[], const IndexType n );

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
         *  @param[in,out]array contains  values and later the differences
         *  @param[in]    n is the number of values, array must contain at least this number of vales
         *  @returns      0
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
        /** Stable sorting of values in array in descending order.
         *
         *  @param[in,out] array are the values to be sorted
         *  @param[in,out] perm, where perm[i] has the value of the original position
         *  @param[in]    n is the number of values to be sorted
         *
         *  \code
         *           array =   1  4   1  8  5  7
         *           perm  =   0  1   2  3  4  5
         +
         *           array =   8  7   5  4  1  1
         *           perm  =   3  5   4  1  0  2
         *  \endcode
         *
         *  Note: if the values to be sorted are int values in a certain range like 0, ..., nb - 1, 
         *        sorting might be more efficient with countBuckets, sortInBuckets
         */

        typedef void ( *FuncType ) ( ValueType array[], IndexType perm[], const IndexType n );

        static const char* getId()
        {
            return "Utils.sort";
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

    template<typename ValueType>
    struct countNonZeros
    {
        /** Count the non-zero elements in an array, used to allocate data for sparse version.
         *
         *  @param[in] denseArray are the values
         *  @param[in] eps        threshold when a value is to be considered as non-zero
         *  @param[in] n          number of elements in the dense array
         *  @returns   number of non-zero elements in denseArray
         */

        typedef IndexType ( *FuncType ) ( const ValueType denseArray[], const IndexType n, const ValueType eps );

        static const char* getId()
        {
            return "Utils.countNonZeros";
        }
    };

    template<typename ValueType>
    struct compress
    {
        /** Build sparse array and sparse indexes from dense array
         *
         *  @param[out] sparseArray     array with non-zero values
         *  @param[out] sparseIndexes  indexes of the non-zero values of input array
         *  @param[in]  denseArray      array with dense values
         *  @param[in]  n               number of elements in the dense array
         *  @returns    number of non-zero elements in denseArray
         *
         *  Note: sparseArray and sparseIndexes must have been allocated with the  correct size before
         */

        typedef IndexType ( *FuncType ) (
            ValueType sparseArray[],
            IndexType sparseIndexes[],
            const ValueType denseArray[],
            const IndexType n,
            const ValueType eps );

        static const char* getId()
        {
            return "Utils.compress";
        }
    };

    template<typename ValueType>
    struct pow
    {
        /** @brief Calculates elementwise pow(array1, array2)
         *
         *  @param[out] result is the result array
         *  @param[in]  array1 is the base array
         *  @param[in]  array2 is the exponent array
         *  @param      n      is the number of entries
         */

        typedef void ( *FuncType ) ( ValueType array1[], const ValueType array2[], const IndexType n );

        static const char* getId()
        {
            return "Util.pow";
        }
    };

    template<typename ValueType>
    struct powBase
    {
        /** @brief Calculates elementwise pow(array1, array2)
         *
         *  @param[out] result is the result array
         *  @param[in]  array is the exponent array
         *  @param[in]  base  is the base
         *  @param      n     is the number of entries
         */

        typedef void ( *FuncType ) ( ValueType array[], ValueType base, const IndexType n );

        static const char* getId()
        {
            return "Util.powBase";
        }
    };

    template<typename ValueType>
    struct powExp
    {
        /** @brief Calculates elementwise pow(array1, array2)
         *
         *  @param[out] result is the result array
         *  @param[in]  array is the base array
         *  @param[in]  exp   is the exponent
         *  @param      n     is the number of entries
         */

        typedef void ( *FuncType ) ( ValueType array[], ValueType exp, const IndexType n );

        static const char* getId()
        {
            return "Util.powExp";
        }
    };
};

} /* end namespace utilskernel */

} /* end namespace scai */

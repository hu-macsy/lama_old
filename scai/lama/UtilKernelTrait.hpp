/**
 * @file UtilKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA utilities provided as kernels.
 * @author Thomas Brandes
 * @date 03.04.2013
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
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
        static const char* getId() { return "Util.validIndexes"; }
    };

    /** @brief Trait for register kernel function sum that sums elements of an array
     *
     *  @tparam ValueType specifies the value type used in the sum reduction.
     */
    template <typename ValueType>
    struct sum
    {
        /** @brief Sum n contiguously stored values.
         *
         *  @param[in] array is an array of values
         *  @param[in] n is the size of array
         *  @return sum of all values in array
         */
        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId () { return "Util.sum"; }
    };

    template <typename ValueType>
    struct maxval
    {
        /** @brief Find maximal value of n contiguously stored values.
         *
         *  @param[in] array is an array of values
         *  @param[in] n is the size of array
         *  @return maximum of all values in array
         */

        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId() { return "Util.maxval"; }
    };

    template <typename ValueType>
    struct absMaxVal
    {
        /** @brief Find absolute maximal value of n contiguously stored values. */

        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId() { return "Util.absMaxVal"; }
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
        static const char* getId() { return "Util.absMaxDiffVal"; }
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
        static const char* getId() { return "Util.isSorted"; }
    };

    /** @brief Structure with functiooń pointer type defintions for setter methods.
     *
     *  @tparam ValueType specifies the value type used in the set operations.
     */

    template<typename ValueType>
    struct setVal
    {
        /** Set all elements of a contiguous array with a value. */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n, const ValueType val );
        static const char* getId() { return "Util.setVal"; }
    };

    template<typename ValueType>
    struct setOrder
    {
        /** Set all elements of a contiguous array with its order number 0, 1, 2, ... */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n );
        static const char* getId() { return "Util.setOrder"; }
    };

    template<typename ValueType>
    struct getValue
    {
        typedef ValueType ( *FuncType ) ( const ValueType* array, const IndexType i );
        static const char* getId() { return "Util.getValue"; }
    };

    template<typename ValueType1, typename ValueType2>
    struct set
    {
        /** Set out[i] = in[i],  0 <= i < n */

        typedef void ( *FuncType ) ( ValueType1 out[], const ValueType2 in[], const IndexType n );
        static const char* getId() { return "Util.set"; }
    };

    template<typename ValueType1, typename ValueType2>
    struct setScale
    {
        /** @brief scaled array assignment, out = in * value
         *
         *  Set out[i] = scale * in[i],  0 <= i < n
         *
         *  @param[in,out]  outValues  is the output array
         *  @param[in]      scaleValue scaling factor
         *  @param[in,out]  inValues   is the array with entries to scale
         *  @param[in]      n          is the number of entries
         */
        typedef void ( *FuncType ) ( ValueType1 outValues[],
                        const ValueType1 scaleValue,
                        const ValueType2 inValues[],
                        const IndexType n );
 
         static const char* getId() { return "Util.setScale"; } 
    };

    template<typename ValueType1, typename ValueType2>
    struct setGather
    {
        /** Set out[i] = in[ indexes[i] ],  \f$0 \le i < n\f$ */

        typedef void ( *FuncType ) ( ValueType1 out[],
                        const ValueType2 in[],
                        const IndexType indexes[],
                        const IndexType n );

        static const char* getId() { return "Util.setGather"; } 
    };

    template<typename ValueType1, typename ValueType2>
    struct setScatter
    {
        /** Set out[ indexes[i] ] = in [i] */

        typedef void ( *FuncType ) ( ValueType1 out[],
                        const IndexType indexes[],
                        const ValueType2 in[],
                        const IndexType n );

        static const char* getId() { return "Util.setScatter"; }
    };

    template<typename ValueType>
    struct invert
    {
        /** @brief Set array[i] = 1.0 / array[i],  0 <= i < n
         *
         *  @param[in,out] array is the array to invert
         *  @param         n     is the number of entries to invert
         */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n );

        static const char* getId() { return "Util.invert"; }
    };

    template<typename ValueType>
    struct scale
    {
        /** @brief scale array of values with a value in place
         *
         *  @param[in,out]  values is the array with entries to scale
         *  @param[in]      value  is the scaling factor
         *  @param[in]      n      is the number of entries in values
         */
        typedef void ( *FuncType ) ( ValueType values[],
                        const ValueType value,
                        const IndexType n );

        static const char* getId() { return "Util.scale"; }
    };
};

} /* end namespace lama */

} /* end namespace scai */

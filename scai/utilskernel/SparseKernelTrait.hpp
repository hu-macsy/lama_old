/**
 * @file SparseKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA utilities involving sparse arrays
 * @author Thomas Brandes
 * @date 26.01.2017
 */
#pragma once

// for dll_import

#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/utilskernel/UnaryOp.hpp>

namespace scai
{

/** Namespace for utilities on heterogeneous arrays (HArray) and derived class LArray */

namespace utilskernel
{

/** Structure just to group traits for all Utils kernels.
 *
 *  This struct does not contain any data at all.
 *  Therefore it could also be a namespace but it is more convenient as
 *  each trait must always be used qualified: SparseKernelTrait::utiliy
 */

struct SparseKernelTrait
{
    template<typename ValueType>
    struct countNonZeros
    {
        /** Count the non-zero elements in an array, used to allocate data for sparse version.
         *
         *  @param[in] denseArray are the values
         *  @param[in] n          number of elements in the dense array
         *  @param[in] eps        threshold when a value is to be considered as non-zero
         *  @returns   number of non-zero elements in denseArray
         */

        typedef IndexType ( *FuncType ) ( const ValueType denseArray[], const IndexType n, const ValueType eps );

        static const char* getId()
        {
            return "Utils.countNonZeros";
        }
    };

    template<typename TargetValueType, typename SourceValueType>
    struct compress
    {
        /** Build sparse array and sparse indexes from dense array
         *
         *  @param[out] sparseArray     array with non-zero values
         *  @param[out] sparseIndexes   indexes of the non-zero values of input array
         *  @param[in]  denseArray      array with dense values
         *  @param[in]  n               number of elements in the dense array
         *  @param[in]  eps             threshold when a value is still considered as zero
         *  @returns    number of non-zero elements in denseArray
         *
         *  Note: the returned value is exactly the same as countNonZeros( denseArray, n, eps )
         *  Note: sparseArray and sparseIndexes must have been allocated with the correct size before
         */

        typedef IndexType ( *FuncType ) (
            TargetValueType sparseArray[],
            IndexType sparseIndexes[],
            const SourceValueType denseArray[],
            const IndexType n,
            const SourceValueType eps );

        static const char* getId()
        {
            return "Utils.compress";
        }
    };

    struct countAddSparse
    {
        /** Count number of indexes for the union of two sparse index sets.
         *
         *  @param[in]  indexes1   non-zero indexes first array
         *  @param[in]  n1         number of non-zero indexes first arra1
         *  @param[in]  indexes2   non-zero indexes first array
         *  @param[in]  n2         number of non-zero indexes second array
         *  @return number of indexes for the union of indexes1 and indexes2
         *
         *  The returned value can be used to allocate the correct size for the sparse array when
         *  two other sparse arrays are added.
         *
         *  It is assumed that the values in indexes1 and indexes2 are sorted in increasing order.
         */

        typedef IndexType ( *FuncType ) (
            const IndexType indexes1[],
            const IndexType n1,
            const IndexType indexes2[],
            const IndexType n2 );

        static const char* getId()
        {
            return "Utils.countAddSparse";
        }
    };

    template<typename ValueType>
    struct addSparse
    {
        /** Add two sparse arrays
         *
         *  @param[out] indexes    non-zero indexes result array
         *  @param[out] values     non-zero values result array
         *  @param[in]  indexes1   non-zero indexes first array
         *  @param[in]  values1    non-zero values first array
         *  @param[in]  n1         number of non-zero indexes first arra1
         *  @param[in]  indexes2   non-zero indexes first array
         *  @param[in]  values2    non-zero values second array
         *  @param[in]  n2         number of non-zero indexes second array
         *  @return     number of non-zero indexes in result array
         *
         *  The returned value must be exactly the same as countAddSparse( indexes1, n1, indexes2, n2 ).
         *  The arrays indexes and values must have been allocated at least with this size.
         */

        typedef IndexType ( *FuncType ) (
            IndexType indexes[],
            ValueType values[],
            const IndexType indexes1[],
            const ValueType values1[],
            const IndexType n1,
            const IndexType indexes2[],
            const ValueType values2[],
            const IndexType n2 );

        static const char* getId()
        {
            return "Utils.addSparse";
        }
    };
};

} /* end namespace utilskernel */

} /* end namespace scai */

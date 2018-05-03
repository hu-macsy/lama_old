/**
 * @file SectionKernelTrait.hpp
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
 * @brief Struct with traits for LAMA utilities to operate on multidimensional array sections
 * @author Thomas Brandes
 * @date 03.04.2013
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>

namespace scai
{

namespace utilskernel
{

/** Structure just to group traits for all section kernels.
 *
 *  This struct does not contain any data at all.
 *  Therefore it could also be a namespace but it is more convenient as
 *  each trait must always be used qualified: SectionKernelTrait::utiliy
 */

struct SectionKernelTrait
{
    /** @brief Trait for register kernel function reduce that reduces elements of an array
     *
     *  @tparam ValueType specifies the value type used in the reduction.
     */
    template <typename ValueType>
    struct assign
    {
        /** @brief assign one section to another
         *
         *  \code
         *     targetSection <op>= sourceSection;
         *  \endcode
         */
        typedef void ( *FuncType ) ( ValueType targetSection[],
                                     const IndexType nDims,
                                     const IndexType sizes[],
                                     const IndexType targetDistances[],
                                     const ValueType sourceSection[],
                                     const IndexType sourceDistances[],
                                     const common::BinaryOp op,
                                     const bool swapOperands );
        static const char* getId()
        {
            return "Section.assign";
        }
    };

    /** @brief Trait for conversion kernel
     *
     *  @tparam SourceValueType specifies the value type used in the reduction.
     *  @tparam TargetValueType specifies the value type used in the reduction.
     */
    template <typename TargetValueType, typename SourceValueType>
    struct unaryOp
    {
        /** @brief assign one section to another
         *
         *  \code
         *     targetSection = unaryOp( sourceSection );
         *  \endcode
         */
        typedef void ( *FuncType ) ( TargetValueType targetSection[],
                                     const IndexType nDims,
                                     const IndexType sizes[],
                                     const IndexType targetDistances[],
                                     const SourceValueType sourceSection[],
                                     const IndexType sourceDistances[],
                                     const common::UnaryOp op );

        static const char* getId()
        {
            return "Section.unaryOp";
        }
    };

    template <typename ValueType>
    struct assignScalar
    {
        /** @brief apply binary operation with a scalar on a sectin
         *
         *
         *  \code
         *     section *= 2.5;
         *     section += 1.2;
         *     section = 0;
         *  \endcode
         */

        typedef void ( *FuncType ) ( ValueType section[],
                                     const IndexType nDims,
                                     const IndexType sizes[],
                                     const IndexType distances[],
                                     ValueType val,
                                     const common::BinaryOp op,
                                     const bool swapOperands );
        static const char* getId()
        {
            return "Section.assign";
        }
    };

    /** @brief Apply unarary operator on a section.
     *
     *  @tparam ValueType specifies the value type of the elements in the section.
     */
    template <typename ValueType>
    struct UnaryOp
    {
        /** @brief assign one section to another
         *
         *  \code
         *     section = abs( section )
         *     section = sin( section )
         *  \endcode
         */

        typedef void ( *FuncType ) ( ValueType section[],
                                     const IndexType nDims,
                                     const IndexType sizes[],
                                     const IndexType distances[],
                                     const common::UnaryOp op );
        static const char* getId()
        {
            return "Section.assign";
        }
    };
};

} /* end namespace utilskernel */

} /* end namespace scai */

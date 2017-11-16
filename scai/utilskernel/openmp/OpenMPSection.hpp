/**
 * @file OpenMPSection.hpp
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
 * @brief Implementation of kernel routines for sections with OpenMP
 * @author Thomas Brandes
 * @date 15.05.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** Copy and operations on sections implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT OpenMPSection
{
public:

    /** OpenMP implementation on host for SectionKernelTrait::assign */

    template<typename ValueType>
    static void assign( 
        ValueType targetSection[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType targetDistances[],
        const ValueType sourceSection[],
        const IndexType sourceDistances[],
        const common::BinaryOp op,
        const bool swapOperands );

    /** OpenMP implementation on host for SectionKernelTrait::assignScalar */

    template <typename ValueType>
    static void assignScalar( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType distances[],
        ValueType val, 
        const common::BinaryOp op,
        const bool swapOperands );

    /** OpenMP implementation on host for SectionKernelTrait::unaryOp */

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp( 
        TargetValueType targetSection[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType targetDistances[],
        const SourceValueType sourceSection[],
        const IndexType sourceDistances[],
        const common::UnaryOp op );

    /** OpenMP implementation on host for SectionKernelTrait::UnaryOp */

    template <typename ValueType>
    static void UnaryOp( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType distances[],
        const common::UnaryOp op );

private:

    /** optimized version of assign with op == COPY */

    template <typename ValueType>
    static void copy2(
        ValueType targetSection[],
        const IndexType sizes[],
        const IndexType targetDistances[],
        const ValueType sourceSection[],
        const IndexType sourceDistances[] );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    struct ArrayKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    template<typename ValueType, typename OtherValueType>
    struct BinOpKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPSection();

    /** Destructor for unregistration. */

    ~OpenMPSection();

    static OpenMPSection guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

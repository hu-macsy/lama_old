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
 * @brief Implementation of general utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
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

/** General utilities of the LAMA Interface implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT OpenMPSection
{
public:

    /** OpenMP implementation on host for SectionKernelTrait::assign */

    template<typename ValueType>
    static void assign( 
        ValueType targetSection[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType targetDifferences[],
        const ValueType sourceSection[],
        const IndexType sourceDifferences[],
        const common::binary::BinaryOp op,
        const bool swapOperands );

    /** OpenMP implementation on host for SectionKernelTrait::assignScalar */

    template <typename ValueType>
    static void assignScalar( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType differences[],
        ValueType val, 
        const common::binary::BinaryOp op,
        const bool swapOperands );

    /** OpenMP implementation on host for SectionKernelTrait::unary */

    template <typename ValueType>
    static void unary( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType differences[],
        const common::unary::UnaryOp op );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    struct ArrayKernels
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

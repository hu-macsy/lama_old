/**
 * @file utilskernel/cuda/CUDASortUtils.hpp
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
 * @brief Implementation of general utilities with CUDA
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

/** General utilities of the LAMA Interface implemented in CUDA  */

class COMMON_DLL_IMPORTEXPORT CUDASortUtils
{
public:

    /** CUDA implementation for UtilKernelTrait::sort */

    template<typename ValueType>
    static void sort(
        IndexType perm[],
        ValueType outValues[],
        const ValueType inValues[],
        const IndexType n,
        const bool ascending );

    /** CUDA implementation for UtilKernelTrait::sortInPlace */

    template<typename ValueType>
    static void sortInPlace(
        IndexType indexes[],
        ValueType values[],
        const IndexType n,
        const bool ascending );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    static void sortPerm( IndexType perm[], const ValueType array[], const IndexType n, bool ascending );

    template<typename ValueType>
    static void sortValues( ValueType array[], const IndexType n, bool ascending );

    template<typename ValueType>
    static void sortBoth( ValueType array[], IndexType perm[], const IndexType n, bool ascending );

    /** Routine that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegArrayKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDASortUtils();

    /** Destructor for unregistration. */

    ~CUDASortUtils();

    /** Static variable for registration at static initialization. */

    static CUDASortUtils guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

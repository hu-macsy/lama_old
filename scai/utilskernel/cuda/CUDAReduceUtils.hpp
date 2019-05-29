/**
 * @file utilskernel/cuda/CUDAReduceUtils.hpp
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
#include <scai/common/CompareOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in CUDA  */

class COMMON_DLL_IMPORTEXPORT CUDAReduceUtils
{
public:

    /** CUDA implementation of UtilKernelTrait::reduce  */

    template<typename ValueType>
    static ValueType reduce(
        const ValueType array[],
        const IndexType n,
        const ValueType zero,
        const common::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::reduce2 */

    template<typename ValueType>
    static ValueType reduce2(
        const ValueType array1[],
        const ValueType array2[],
        const IndexType n,
        const common::BinaryOp binOp,
        const ValueType zero,
        const common::BinaryOp redOp );

    /** CUDA implementation for UtilKernelTrait::allCompare */

    template<typename ValueType>
    static bool allCompare(
        const ValueType array1[],
        const ValueType array2[],
        const IndexType n,
        const common::CompareOp op );

    /** CUDA implementation for UtilKernelTrait::allCompareScalar */

    template<typename ValueType>
    static bool allCompareScalar(
        const ValueType array[],
        const ValueType scalar,
        const IndexType n,
        const common::CompareOp op );

    /** CUDA implementation for UtilKernelTrait::scan */

    template<typename ValueType>
    static ValueType scan(
        ValueType array[],
        const IndexType n,
        const ValueType first,
        const bool exclusive,
        const bool append );

    /** CUDA implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, const common::CompareOp op );

    /** CUDA implementation of UtilKernelTrait::validIndexes, contains bool reduction  */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    static ValueType reduceSum( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceMaxVal( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceMinVal( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceAbsMaxVal( const ValueType array[], const IndexType n, const ValueType zero );

    /** Routine that registers all methods at the kernel registry. */

    struct Registrator
    {   
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    template<typename ValueType>
    struct RegArrayKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDAReduceUtils();

    /** Destructor for unregistration. */

    ~CUDAReduceUtils();

    /** Static variable for registration at static initialization. */

    static CUDAReduceUtils guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

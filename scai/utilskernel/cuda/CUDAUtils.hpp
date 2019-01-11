/**
 * @file utilskernel/cuda/CUDAUtils.hpp
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

class COMMON_DLL_IMPORTEXPORT CUDAUtils
{
public:

    /** CUDA implementation of UtilKernelTrait::setVal  */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const common::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::scaleVectorAddScalar */

    template<typename ValueType>
    static void scaleVectorAddScalar( ValueType array1[], const ValueType array2[], const IndexType n, const ValueType alpha, const ValueType beta );

    /** CUDA implementation of UtilKernelTrait::getValue  */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** CUDA implementation for UtilKernelTrait::unaryOp */

    template<typename ValueType>
    static void unaryOp( ValueType out[], const ValueType in[], const IndexType n, const common::UnaryOp op );

    /** CUDA implementation for UtilKernelTrait::binaryOp */

    template<typename ValueType>
    static void binaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const common::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::binaryOpScalar */

    template<typename ValueType>
    static void binaryOpScalar(
        ValueType out[],
        const ValueType in[],
        const ValueType value,
        const IndexType n,
        const common::BinaryOp op,
        const bool swapScalar );

    /** CUDA implementation for UtilKernelTrait::scatterVal */

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    template<typename ValueType, typename OtherValueType>
    struct RegistratorVO
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDAUtils();

    /** Destructor for unregistration. */

    ~CUDAUtils();

    /** Static variable for registration at static initialization. */

    static CUDAUtils guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

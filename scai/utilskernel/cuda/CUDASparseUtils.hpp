/**
 * @file utilskernel/cuda/CUDASparseUtils.hpp
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

class COMMON_DLL_IMPORTEXPORT CUDASparseUtils
{
public:

    /** CUDA implementation of UtilKernelTrait::setOrder  */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::Setter::setSequence */

    template<typename ValueType>
    static void setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::set */

    template<typename ValueType, typename otherValueType>
    static void set( ValueType out[], const otherValueType in[], const IndexType n, const common::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::setSection */

    template<typename ValueType, typename otherValueType>
    static void setSection( ValueType out[], const IndexType inc_out,
                            const otherValueType in[], const IndexType inc_in, const IndexType n, const common::BinaryOp op );


    /** CUDA implementation for UtilKernelTrait::setGather */

    template<typename ValueType, typename otherValueType>
    static void setGather(
        ValueType out[],
        const otherValueType in[],
        const IndexType indexes[],
        const common::BinaryOp op,
        const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScatter */

    template<typename ValueType, typename otherValueType>
    static void setScatter(
        ValueType out[],
        const IndexType indexes[],
        const bool unique,
        const otherValueType in[],
        const common::BinaryOp op,
        const IndexType n );

    /** CUDA implementation of UtilKernelTrait::countNonZeros */

    template<typename ValueType>
    static IndexType countNonZeros( 
        const ValueType denseArray[], 
        const IndexType n, 
        const ValueType zero, 
        ValueType eps );

    /** CUDA implementation of UtilKernelTrait::compress */

    template<typename TargetValueType, typename SourceValueType>
    static IndexType compress(
        TargetValueType sparseArray[],
        IndexType sparseIndexes[],
        const SourceValueType denseArray[],
        const IndexType n,
        const SourceValueType zero,
        const SourceValueType eps );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    CUDASparseUtils();

    /** Destructor for unregistration. */

    ~CUDASparseUtils();

    /** Static variable for registration at static initialization. */

    static CUDASparseUtils guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

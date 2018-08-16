/**
 * @file CUDASection.hpp
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
 * @brief Implementation of kernel routines for sections with CUDA
 * @author Thomas Brandes
 * @date 17.05.2017
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

class COMMON_DLL_IMPORTEXPORT CUDASection
{
public:

    /** CUDA implementation on host for SectionKernelTrait::assign */

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

    /** CUDA implementation on host for SectionKernelTrait::assignScalar */

    template <typename ValueType>
    static void assignScalar( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType distances[],
        ValueType val, 
        const common::BinaryOp op,
        const bool swapOperands );

    /** CUDA implementation on host for SectionKernelTrait::unaryOp */

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp( 
        TargetValueType targetSection[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType targetDistances[],
        const SourceValueType sourceSection[],
        const IndexType sourceDistances[],
        const common::UnaryOp op );

    /** CUDA implementation on host for SectionKernelTrait::UnaryOp */

    template <typename ValueType>
    static void UnaryOp( 
        ValueType section[],
        const IndexType nDims,
        const IndexType sizes[],
        const IndexType distances[],
        const common::UnaryOp op );

private:

    template<typename ValueType>
    static void assign0(
        ValueType targetSection[],
        const ValueType sourceSection[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assign1(
        ValueType targetSection[],
        const ValueType sourceSection[],
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assign2(
        ValueType targetSection[],
        const ValueType sourceSection[],
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assign3(
        ValueType targetSection[],
        const ValueType sourceSection[],
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assign4(
        ValueType targetSection[],
        const ValueType sourceSection[],
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assignScalar0(
        ValueType targetSection[],
        const ValueType val,
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assignScalar1(
        ValueType targetSection[],
        const ValueType val,
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assignScalar2(
        ValueType targetSection[],
        const ValueType val,
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assignScalar3(
        ValueType targetSection[],
        const ValueType val,
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename ValueType>
    static void assignScalar4(
        ValueType targetSection[],
        const ValueType val,
        const IndexType sizes[],
        const common::BinaryOp op,
        const bool swapOperands );

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp1(
        TargetValueType targetSection[],
        const SourceValueType sourceSection[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp2(
        TargetValueType targetSection[],
        const SourceValueType sourceSection[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp3(
        TargetValueType targetSection[],
        const SourceValueType sourceSection[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename TargetValueType, typename SourceValueType>
    static void unaryOp4(
        TargetValueType targetSection[],
        const SourceValueType sourceSection[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename ValueType>
    static void UnaryOp1(
        ValueType section[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename ValueType>
    static void UnaryOp2(
        ValueType section[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename ValueType>
    static void UnaryOp3(
        ValueType section[],
        const IndexType sizes[],
        const common::UnaryOp op );

    template<typename ValueType>
    static void UnaryOp4(
        ValueType section[],
        const IndexType sizes[],
        const common::UnaryOp op );

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

    CUDASection();

    /** Destructor for unregistration. */

    ~CUDASection();

    static CUDASection guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

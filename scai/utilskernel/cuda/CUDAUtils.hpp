/**
 * @file utilskernel/cuda/CUDAUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of general utilities with CUDA
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/utilskernel/UnaryOp.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in CUDA  */

class COMMON_DLL_IMPORTEXPORT CUDAUtils
{
public:

    /** CUDA implementation of UtilKernelTrait::validIndexes  */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** CUDA implementation of UtilKernelTrait::reduce  */

    template<typename ValueType>
    static ValueType reduce(
        const ValueType array[],
        const IndexType n,
        const ValueType zero,
        const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::reduce2 */

    template<typename ValueType>
    static ValueType reduce2(
        const ValueType array1[],
        const ValueType array2[],
        const IndexType n,
        const binary::BinaryOp binOp,
        const ValueType zero,
        const binary::BinaryOp redOp );

    /** CUDA implementation of UtilKernelTrait::setVal  */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::scaleVectorAddScalar */

    template<typename ValueType>
    static void scaleVectorAddScalar( ValueType array1[], const ValueType array2[], const IndexType n, const ValueType alpha, const ValueType beta );

    /** CUDA implementation of UtilKernelTrait::setOrder  */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::Setter::setSequence */

    template<typename ValueType>
    static void setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );

    /** CUDA implementation of UtilKernelTrait::getValue  */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** CUDA implementation UtilKernelTrait::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    /** CUDA implementation for UtilKernelTrait::set */

    template<typename ValueType, typename otherValueType>
    static void set( ValueType out[], const otherValueType in[], const IndexType n, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::setSection */

    template<typename ValueType, typename otherValueType>
    static void setSection( ValueType out[], const IndexType inc_out, 
                            const otherValueType in[], const IndexType inc_in, const IndexType n, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::unaryOp */

    template<typename ValueType>
    static void unaryOp( ValueType out[], const ValueType in[], const IndexType n, const unary::UnaryOp op );

    /** CUDA implementation for UtilKernelTrait::binaryOp */

    template<typename ValueType>
    static void binaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::binaryOpScalar1 */

    template<typename ValueType>
    static void binaryOpScalar1( ValueType out[], const ValueType value, const ValueType in[], const IndexType n, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::binaryOpScalar2 */

    template<typename ValueType>
    static void binaryOpScalar2( ValueType out[], const ValueType in[], const ValueType value, const IndexType n, const binary::BinaryOp op );

    /** CUDA implementation for UtilKernelTrait::setGather */

    template<typename ValueType, typename otherValueType>
    static void setGather( 
        ValueType out[], 
        const otherValueType in[], 
        const IndexType indexes[], 
        const utilskernel::binary::BinaryOp op,
        const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScatter */

    template<typename ValueType, typename otherValueType>
    static void setScatter( 
        ValueType out[], 
        const IndexType indexes[], 
        const otherValueType in[], 
        const binary::BinaryOp op, 
        const IndexType n );

    /** CUDA implementation for UtilKernelTrait::scatterVal */

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::scan */

    template<typename ValueType>
    static ValueType scan( ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::sort */

    template<typename ValueType>
    static void sort( ValueType array[], IndexType perm[], const IndexType n, const bool ascending );

    /** CUDA implementation for UtilKernelTrait::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** CUDA implementation of UtilKernelTrait::countNonZeros */

    template<typename ValueType>
    static IndexType countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType eps );

    /** CUDA implementation of UtilKernelTrait::compress */

    template<typename ValueType>
    static IndexType compress(
        ValueType sparseArray[],
        IndexType sparseIndexes[],
        const ValueType denseArray[],
        const IndexType n,
        const ValueType eps );

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

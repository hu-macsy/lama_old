/**
 * @file CUDAUtils.hpp
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

#include <scai/utilskernel/ReductionOp.hpp>
#include <scai/utilskernel/ElementwiseOp.hpp>

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

    /** CUDA implementation for UtilKernelTrait::vectorScale */

    template<typename ValueType>
    static void vectorScale( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n );

    /*  CUDA implementation of UtilKernelTrait::validIndexes  */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /*  CUDA implementation of UtilKernelTrait::reduce  */

    template<typename ValueType>
    static ValueType reduce( const ValueType array[], const IndexType n, const reduction::ReductionOp op );

    /*  CUDA implementation of UtilKernelTrait::setVal  */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const reduction::ReductionOp op );

    /*  CUDA implementation of UtilKernelTrait::setOrder  */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::Setter::setSequence */

    template<typename ValueType>
    static void setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );

    /*  CUDA implementation of UtilKernelTrait::getValue  */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** CUDA implementation for UtilKernelTrait::scale. */

    template<typename ValueType>
    static void scale( ValueType values[], const ValueType value, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScale. */

    template<typename ValueType, typename otherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const otherValueType inValues[],
        const IndexType n );

    /** CUDA function implements UtilKernelTrait::absMaxVal */

    template<typename ValueType>
    static ValueType absMaxVal( const ValueType array[], const IndexType n );

    /** CUDA function implements UtilKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::copysign */

    template<typename ValueType>
    static void copysign( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    /** CUDA implementation for UtilKernelTrait::set */

    template<typename ValueType, typename otherValueType>
    static void set( ValueType out[], const otherValueType in[], const IndexType n, const reduction::ReductionOp op );

   /** CUDA implementation for UtilKernelTrait::execElementwise */

    template<typename ValueType>
    static void execElementwise( ValueType array[], const IndexType n, const elementwise::ElementwiseOp op );

   /** CUDA implementation for UtilKernelTrait::pow */

    template<typename ValueType>
    static void pow( ValueType array1[], const ValueType array2[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::powBasw */

    template<typename ValueType>
    static void powBase( ValueType array[], const ValueType base, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::powExp */

    template<typename ValueType>
    static void powExp( ValueType array[], const ValueType exp, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setGather, out[i]] = in[ indexes[i] ] */

    template<typename ValueType, typename otherValueType>
    static void setGather( ValueType out[], const otherValueType in[], const IndexType indexes[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::setScatter, out[ indexes[i] ] op= in [i] */

    template<typename ValueType, typename otherValueType>
    static void setScatter( ValueType out[], const IndexType indexes[], const otherValueType in[], const reduction::ReductionOp op, const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::scatterVal */

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

    /** CUDA implementation for UtilKernelTrait::scan */

    template<typename ValueType>
    static ValueType scan( ValueType array[], const IndexType n );

    /** CUDA implementation for UtilKernelTrait::sort */

    template<typename ValueType>
    static void sort( ValueType array[], IndexType perm[], const IndexType n );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    static ValueType reduceSum( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMaxVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMinVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceAbsMaxVal( const ValueType array[], const IndexType n );

    /** Routine that registers all methods at the kernel registry. */

    struct Registrator
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag ); 
    };

    template<typename ValueType>
    struct RegNumericKernels
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

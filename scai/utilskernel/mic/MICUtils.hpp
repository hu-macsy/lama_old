/**
 * @file MICUtils.hpp
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
 * @brief Implementation of general utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/utilskernel/UnaryOp.hpp>

// logging
#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in MIC  */

class COMMON_DLL_IMPORTEXPORT MICUtils
{
public:

    /*  This method is an implementation of UtilKernelTrait::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** MIC implementation for UtilKernelTrait::Reductions::reduce */

    template<typename ValueType>
    static ValueType reduce( const ValueType array[], const IndexType n, const binary::BinaryOp op );

    /** MIC implementation for UtilKernelTrait::Setter::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const binary::BinaryOp op );

    /** MIC implementation for UtilKernelTrait::Setter::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** MIC implementation for UtilKernelTrait::Setter::setSequence */

    template<typename ValueType>
    static void setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );

    /** MIC implementation for UtilKernelTrait::getValue */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** MIC implementation for UtilKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** MIC implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    template<typename ValueType1, typename ValueType2>
    static void set( ValueType1 out[], const ValueType2 in[], const IndexType n, const binary::BinaryOp op );

    /** MIC implementation for UtilKernelTrait::setGather, out[i]] op = in[ indexes[i] ] */

    template<typename ValueType, typename otherValueType>
    static void setGather(
        ValueType out[],
        const otherValueType in[],
        const IndexType indexes[],
        const utilskernel::binary::BinaryOp op,
        const IndexType n );

    /** Set out[ indexes[i] ] = in [i] */

    template<typename ValueType1, typename ValueType2>
    static void setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const binary::BinaryOp op, const IndexType n );

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

    /** MIC implementation for UtilKernelTrait::applyUnaryOp */

    template<typename ValueType>
    static void applyUnaryOp( ValueType out[], const ValueType in[], const IndexType n, const unary::UnaryOp op );

    /** MIC implementation for UtilKernelTrait::applyBinaryOp */

    template<typename ValueType>
    static void applyBinaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const binary::BinaryOp op );

    /** MIC implementation for UtilKernelTrait::applyBinaryOpScalar1 */

    template<typename ValueType>
    static void applyBinaryOpScalar1( ValueType out[], const ValueType value, const ValueType in[], const IndexType n, const binary::BinaryOp op );

    /** MIC implementation for UtilKernelTrait::applyBinaryOpScalar2 */

    template<typename ValueType>
    static void applyBinaryOpScalar2( ValueType out[], const ValueType in[], const ValueType value, const IndexType n, const binary::BinaryOp op );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    template<typename ValueType>
    static ValueType reduceSum( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMinVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMaxVal( const ValueType array[], const IndexType n );

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

    /** Helper class for (un) registration of kernel routines at static initialization. */

    class RegisterGuard
    {
    public:
        RegisterGuard();
        ~RegisterGuard();
    };

    static RegisterGuard guard;  // registration of kernels @ static initialization

};

} /* end namespace utilskernel */

} /* end namespace scai */

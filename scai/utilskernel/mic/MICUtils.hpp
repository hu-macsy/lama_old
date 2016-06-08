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
#include <scai/utilskernel/ReductionOp.hpp>

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

    /** MIC implementation for UtilKernelTrait::Copy::setScale */

    template<typename ValueType, typename OtherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const OtherValueType inValues[],
        const IndexType n );

    /*  This method is an implementation of UtilKernelTrait::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** MIC implementation for UtilKernelTrait::Reductions::reduce */

    template<typename ValueType>
    static ValueType reduce( const ValueType array[], const IndexType n, const reduction::ReductionOp op );

    /** MIC implementation for UtilKernelTrait::Setter::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const reduction::ReductionOp op );

    /** MIC implementation for UtilKernelTrait::Setter::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** MIC implementation for UtilKernelTrait::getValue */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** MIC implementation for UtilKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** MIC implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    template<typename ValueType1,typename ValueType2>
    static void set( ValueType1 out[], const ValueType2 in[], const IndexType n, const reduction::ReductionOp op );

    /** Set out[i] = in[ indexes[i] ],  0 <= i < n */

    template<typename ValueType1, typename ValueType2>
    static void setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n );

    /** Set out[ indexes[i] ] = in [i] */

    template<typename ValueType1, typename ValueType2>
    static void setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n );

    /** MIC implementation for UtilKernelTrait::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

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

    SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )


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

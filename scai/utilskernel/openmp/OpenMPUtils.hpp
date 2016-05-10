/**
 * @file OpenMPUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
#include <scai/common/ReductionOp.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT OpenMPUtils
{
public:

    /** OpenMP implementation for UtilKernelTrait::conj */

    template<typename ValueType>
    static void conj( ValueType mValues[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::setScale */

    template<typename ValueType,typename OtherValueType>
    static void setScale(
        ValueType outValues[],
        const ValueType value,
        const OtherValueType inValues[],
        const IndexType n );

    /*  This method is an implementation of UtilKernelTrait::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** OpenMP implementation for UtilKernelTrait::reduce */

    template<typename ValueType>
    static ValueType reduce( const ValueType array[], const IndexType n, const common::reduction::ReductionOp op );

    /** OpenMP implementation for UtilKernelTrait::Setter::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const common::reduction::ReductionOp op );

    /** OpenMP implementation for UtilKernelTrait::Setter::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** OpenMP implementation for UtilKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, bool acending );

    /** OpenMP implementation for UtilKernelTrait::set */

    template<typename ValueType1,typename ValueType2>
    static void set( ValueType1 out[], const ValueType2 in[], const IndexType n, const common::reduction::ReductionOp op );

    /** OpenMP implementation for UtilKernelTrait::setGather */

    template<typename ValueType1,typename ValueType2>
    static void setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::scatterVal */

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::setScatter */

    template<typename ValueType1,typename ValueType2>
    static void setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::invert */

    template<typename ValueType>
    static void invert( ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::scan */

    template<typename ValueType>
    static ValueType scan( ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::sort, uses bucket sort */

    template<typename ValueType>
    static void sort( ValueType array[], IndexType perm[], const IndexType n );

    /** Compute the inverse permutation as specified in UtilKernelTrait::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

private:

    template<typename ValueType>
    static ValueType reduceSum( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMaxVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceMinVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType reduceAbsMaxVal( const ValueType array[], const IndexType n );

    template<typename ValueType>
    static ValueType scanSerial( ValueType array[], const IndexType numValues );

    template<typename ValueType>
    static ValueType scanParallel( PartitionId numThreads, ValueType array[], const IndexType numValues );

    /** OpenMP implementation of UtilsKernelTrait::countNonZeros */

    template<typename ValueType>
    static IndexType countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType eps );

    /** OpenMP implementation of UtilsKernelTrait::compress */

    template<typename ValueType>
    static IndexType compress( 
        ValueType sparseArray[], 
        IndexType sparseIndexes[], 
        const ValueType denseArray[], 
        const IndexType n, 
        const ValueType eps );

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

    /** Constructor for registration. */

    OpenMPUtils();

    /** Destructor for unregistration. */

    ~OpenMPUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace utilskernel */

} /* end namespace scai */

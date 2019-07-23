/**
 * @file CUDADIAUtils.hpp
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
 * @brief Implementation of DIA utilities with CUDA
 * @author Thomas Brandes
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides CUDA parallelized routines needed for DIA format.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDADIAUtils
{
public:

    /** Implementation for DIAKernelTrait::normalGEMV with CUDA on GPUs */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffsets[],
        const ValueType diaValues[],
        const common::MatrixOp op );

    /** Implementation for DIAKernelTrait::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType n,
        const IndexType numDiagonals,
        const IndexType diaOffset[],
        const ValueType diaValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega );

    /** Implementation for DIAKernelTrait::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const ValueType diagonal[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numDiagonals,
        const IndexType diaOffset[],
        const ValueType diaValues[],
        const ValueType oldSolution[],
        const ValueType omega );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDADIAUtils();

    /** Destructor for unregistration. */

    ~CUDADIAUtils();

    /** Static variable for registration at static initialization. */

    static CUDADIAUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */

/**
 * @file CUDACOOUtils.hpp
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
 * @brief Implementation of COO utilities with CUDA
 * @author Thomas Brandes
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides CUDA parallelized routines needed for COO format.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDACOOUtils
{
public:

    /** Implementation for COOKernelTrait::offsets2ia with CUDA on GPUs */

    static void offsets2ia(
        IndexType cooIA[],
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType numRows );

    /** Help routine that computes an offsets array for sorted array of COO row indexes */

    static void ia2offsets(
        IndexType csrIA[],
        const IndexType numRows,
        const IndexType cooIA[],
        const IndexType numValues );

    /** Implementation for COOKernelTrait::normalGEMV with CUDA on GPUs */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const common::MatrixOp op );

    /** Implementation for COOKernelTrait::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType cooNumValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega,
        const IndexType numRows );

    /** Implementation for COOKernelTrait::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType cooNumValues,
        const IndexType cooIA[],
        const IndexType cooJA[],
        const ValueType cooValues[],
        const ValueType localDiagonal[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numRows );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Struct for registration of methods without template arguments */

    struct Registrator
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

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

    CUDACOOUtils();

    /** Destructor for unregistration. */

    ~CUDACOOUtils();

    /** Static variable for registration at static initialization. */

    static CUDACOOUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */

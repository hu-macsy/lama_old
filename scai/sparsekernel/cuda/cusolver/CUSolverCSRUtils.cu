/**
 * @file CUSolverCSRUtils.cu
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
 * @brief Implementation of some CSR routines with CUSolver library
 * @author Thomas Brandes
 * @date 11.06.2013
 */

// hpp
#include <scai/sparsekernel/cuda/cusolver/CUSolverCSRUtils.hpp>

// local library
#include <scai/sparsekernel/cuda/cusolver/CUSOLVERWrapper.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// CUDA
#include <cuda.h>
#include <cusparse_v2.h>

#if ( CUDART_VERSION >= 7050 )

namespace scai
{

using tasking::CUDAStreamSyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUSolverCSRUtils::logger, "CUDA.CUSolverCSRUtils" )

/* --------------------------------------------------------------------------- */
/*                          decomposition (cuSolver)                           */
/* --------------------------------------------------------------------------- */


template<typename ValueType>
void CUSolverCSRUtils::decomposition(
    ValueType* const solution,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const ValueType rhs[],
    const IndexType numRows,
    const IndexType nnz,
    const bool isSymmetic )
{
    SCAI_LOG_DEBUG( logger,
                    "decomposition<" << common::getScalarType<ValueType>() << ", matrix numRows = "
                    << numRows << ", nnz = " << nnz )

    SCAI_REGION( "CUDA.cuSolver.decomposition" )

    typedef CUSOLVERTrait::BLASIndexType BLASIndexType;

    if ( common::TypeTraits<IndexType>::stype != common::TypeTraits<BLASIndexType>::stype )
    {
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    // note: SCAI_CHECK_CUDA_ACCESS not required due to getCurrentCUDACtx
    cusolverSpHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSolverSpHandle();

    cusparseMatDescr_t descrA;
    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrA ), "cusparseCreateMatDescr" )

    if ( isSymmetic )
    {
        if ( scai::common::isComplex( scai::common::TypeTraits<ValueType>::stype ) )
        {
            cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_HERMITIAN );
        }
        else
        {
            cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_SYMMETRIC );
        }
    }
    else
    {
        cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL );
    }

    cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO );

    ValueType tol = 1e-10;
    IndexType reorder = 0;
    IndexType singularity = 0;

    // call for QR decomposition
    CUSOLVERWrapper<ValueType>::csrQR( handle, numRows, nnz, descrA, csrValues, csrIA, csrJA,
                                       rhs, tol, reorder, solution, &singularity );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "decomposition" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUSolverCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_INFO( logger, "register CUSolverCSRUtils CUSolver-routines for CUDA at kernel registry ["
                   << flag << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::decomposition<ValueType> >( decomposition, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUSolverCSRUtils::CUSolverCSRUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUSolverCSRUtils::~CUSolverCSRUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUSolverCSRUtils CUSolverCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */

#endif

/**
 * @file MKLCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with MKL
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/external/MKLCSRUtils.hpp>

// local library
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

#include <scai/sparsekernel/external/MKLCSRTrait.hpp>
#include <scai/sparsekernel/external/MKLCSRWrapper.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/unused.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

// extern
#include <mkl_spblas.h>

using namespace scai::utilskernel;

namespace scai
{

using tasking::TaskSyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MKLCSRUtils::logger, "MKL.CSRUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* nnz */,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "MKL.scsrmv" )
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )
    typedef MKLCSRTrait::BLASIndexType BLASIndexType;
    typedef MKLCSRTrait::BLASTrans BLASTrans;
    typedef MKLCSRTrait::BLASMatrix BLASMatrix;
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( common::TypeTraits<IndexType>::stype
            != common::TypeTraits<BLASIndexType>::stype )
    {
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    if ( syncToken )
    {
        SCAI_UNSUPPORTED( "asynchronous execution not supported yet" )
        // ToDo: workaround required as boost::bind supports only up to 9 arguments
    }

    if ( y != result && beta != scai::common::constants::ZERO )
    {
        OpenMPUtils::set( result, y, numRows, utilskernel::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y
    BLASTrans transa = 'n';
    // General, - triangular, Non-Unit, C for zero-indexing
    BLASMatrix matdescra;
    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';
    // const_cast needed, MKL interface does not support it
    MKLCSRWrapper<ValueType>::csrmv( transa, numRows, numColumns, alpha, matdescra, csrValues,
                                     csrJA, csrIA, csrIA + 1, x, beta, result );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    ValueType cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    // Intel MKL supports csr to csc only for square matrices
    if ( numRows == numColumns )
    {
        SCAI_REGION( "MKL.CSRUtils.convertCSR2CSC" )
        SCAI_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )
        IndexType job[] =
        { 0, 0, 0, 0, 0, 1 };
        MKLCSRWrapper<ValueType>::csr2csc( job, numRows, csrValues, csrJA, csrIA, cscValues, cscJA, cscIA );
    }
    else
    {
        OpenMPCSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, csrIA, csrJA, csrValues, numRows, numColumns,
                                        numValues );
    }
}


/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::context::ContextType ctx = common::context::Host;
    using kregistry::KernelRegistry;
    SCAI_LOG_INFO( logger, "register CSRUtils MKL-routines for Host at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

MKLCSRUtils::MKLCSRUtils()
{
    bool useMKL = true;  // by default: use MKL for CSR utils

    common::Settings::getEnvironment( level, "SCAI_USE_MKL" );

    if ( !useMKL )
    {
        return;
    }

    // replace internal OpenMP kernels  

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call( flag );
}

MKLCSRUtils::~MKLCSRUtils()
{
    bool useMKL = true;

    common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );

    if ( !useMKL )
    {
        return;
    }

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

MKLCSRUtils MKLCSRUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */

/**
 * @file MICMKLCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with MICMKL
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/mic/MICMKLCSRUtils.hpp>

// local library
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/mic/MICMKLCSRWrapper.hpp>
#include <scai/sparsekernel/external/MKLCSRTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/tracing.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// external
#include <mkl_spblas.h>

namespace scai
{

using namespace hmemo;
using tasking::MICSyncToken;
using utilskernel::MICUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MICMKLCSRUtils::logger, "MIC.MKLCSRUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICMKLCSRUtils::normalGEMV(
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
    // SCAI_REGION( "MIC.MKLscsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    typedef MKLCSRTrait::BLASIndexType BLASIndexType;
    typedef MKLCSRTrait::BLASTrans BLASTrans;
    typedef MKLCSRTrait::BLASMatrix BLASMatrix;

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet" )
    }

    if( y != result && beta != 0 )
    {
        MICUtils::set( result, y, numRows, common::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y

    BLASTrans transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    BLASMatrix matdescra;

    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    MICMKLCSRWrapper<ValueType>::csrmv( transa, numRows, numColumns, alpha, matdescra, csrValues, csrJA, csrIA, csrIA + 1, x, beta, result );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICMKLCSRUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::context::ContextType ctx = common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register CSRUtils MKL-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICMKLCSRUtils::MICMKLCSRUtils()
{
    bool useMKL = true;

   // using MKL for CSR might be disabled explicitly by environment variable

   common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );
   int level = 0;

   if( !useMKL || ( level <= 0 ) )
   {
       return;
   }

   const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

   kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICMKLCSRUtils::~MICMKLCSRUtils()
{
    bool useMKL = true;

   // using MKL for CSR might be disabled explicitly by environment variable

   common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );
   int level = 0;

   if( !useMKL || ( level <= 0 ) )
   {
       return;
   }

   const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

   kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICMKLCSRUtils MICMKLCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */

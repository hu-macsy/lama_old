/**
 * @file MICMKLCSRUtils.cpp
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
#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
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
                   "normalGEMV<" << common::TypeTraits<ValueType>::id() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    typedef MKLCSRTrait::BLASIndexType BLASIndexType;
    typedef MKLCSRTrait::BLASTrans BLASTrans;
    typedef MKLCSRTrait::BLASMatrix BLASMatrix;

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet" )
    }

    if ( y != result && beta != 0 )
    {
        MICUtils::set( result, y, numRows, utilskernel::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y

    BLASTrans transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    const void *iaPtr = csrIA;
    const void *jaPtr = csrJA;
    const void *valPtr = csrValues;
    const void *xPtr = x;
    void *resultPtr = result;
    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( transa, numRows, numColumns, alphaPtr[0:1], valPtr, jaPtr, iaPtr, xPtr, betaPtr[0:1], resultPtr )
    {
        BLASMatrix matdescra;

        ValueType alpha = alphaPtr[0];
        ValueType beta = betaPtr[0];

        matdescra[0] = 'g';
        matdescra[1] = ' ';
        matdescra[2] = 'n';
        matdescra[3] = 'c';

        const IndexType* csrIA = reinterpret_cast<const IndexType*>( iaPtr );
        const IndexType* csrJA = reinterpret_cast<const IndexType*>( jaPtr );
        const ValueType* csrValues = reinterpret_cast<const ValueType*>( valPtr );

        const ValueType* x = reinterpret_cast<const ValueType*>( xPtr );
        ValueType* result = reinterpret_cast<ValueType*>( resultPtr );

        MICMKLCSRWrapper<ValueType>::csrmv( transa, numRows, numColumns, alpha, matdescra, csrValues, csrJA, csrIA, csrIA + 1, x, beta, result );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICMKLCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::context::ContextType ctx = common::context::MIC;

    using kregistry::KernelRegistry;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "] CSRUtils MKL-routines for MIC at kernel registry: "
                            << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICMKLCSRUtils::MICMKLCSRUtils()
{
    bool useMKL = true;   // default is enabled

    // using MKL for CSR might be disabled explicitly by environment variable

    common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );

    if ( !useMKL )
    {
        SCAI_LOG_INFO( logger, "disabled: CSRUtils MKL-routines for MIC at kernel registry" )
        return;
    }

    SCAI_LOG_INFO( logger, "register CSRUtils MKL-routines for MIC at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICMKLCSRUtils::~MICMKLCSRUtils()
{
    bool useMKL = true;

    common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );

    if ( !useMKL )
    {
        SCAI_LOG_INFO( logger, "disabled: CSRUtils MKL-routines for MIC at kernel registry" )
        return;
    }

    SCAI_LOG_INFO( logger, "unregister CSRUtils MKL-routines for MIC at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICMKLCSRUtils MICMKLCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */

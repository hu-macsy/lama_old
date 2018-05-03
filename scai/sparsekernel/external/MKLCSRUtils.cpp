/**
 * @file MKLCSRUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

#include <scai/sparsekernel/external/MKLCSRTrait.hpp>
#include <scai/sparsekernel/external/MKLCSRWrapper.hpp>
#include <scai/sparsekernel/external/PardisoError.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/unused.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

// extern
#include <mkl.h>
#include <mkl_spblas.h>

namespace scai
{

using namespace utilskernel;
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
    const ValueType csrValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "MKL.scsrmv" )
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )
    typedef MKLCSRTrait::BLASIndexType BLASIndexType;
    typedef MKLCSRTrait::BLASTrans BLASTrans;
    typedef MKLCSRTrait::BLASMatrix BLASMatrix;
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( common::TypeTraits<IndexType>::stype != common::TypeTraits<BLASIndexType>::stype )
    {
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    if ( syncToken )
    {
        SCAI_UNSUPPORTED( "asynchronous execution not supported yet" )
        // ToDo: workaround required as boost::bind supports only up to 9 arguments
    }

    if ( y != result && beta != common::Constants::ZERO )
    {
        IndexType n = common::isTranspose( op ) ? numColumns : numRows;

        OpenMPUtils::set( result, y, n, common::BinaryOp::COPY );
    }

    // performs y = alpha * A * x + beta * y

    BLASTrans transa = common::isTranspose( op ) ? 't' : 'n';

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

template<typename ValueType>
void MKLCSRUtils::decomposition(
    ValueType* solution,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const ValueType rhs[],
    const IndexType numRows,
    const IndexType nnz,
    const bool isSymmetric )
{
    SCAI_REGION( "MKL.decomposition" )

    SCAI_LOG_INFO( logger, "decomposition<" << common::TypeTraits<ValueType>::id() << "> of matrix,"
                   << " numRows = " << numRows << ", nnz = " << nnz << ", isSymmetric = " << isSymmetric )

    // dummy variables

    ValueType vDum;
    MKL_INT   iDum;

    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only       */
    /* necessary for the FIRST call of the PARDISO solver.                  */
    /* -------------------------------------------------------------------- */
    void* pt[64];

    for ( int i = 0; i < 64; ++i )
    {
        pt[i] = 0;
    }

    // TODO: check all posibilities of matrix type

    //  1: Real structural symmetrix matrix
    //  2: Real symmetric positive definite (spd) matrix
    // -2: Real symmetric indefinite matrix
    //  3: Complex structural symmetrix matrix
    //  4: Complex hermitian positive definite matrix
    // -4: Complex hermitian indefinite matrix
    //  6: Complex symmetric matrix
    // 11: Real unsymmetric matrix
    // 13: Complex unsymmetric matrix
    MKL_INT mtype;

    if ( common::isComplex( common::TypeTraits<ValueType>::stype ) ) // Complex
    {
        if ( isSymmetric )
        {
            mtype = 6;
            // also may be: 3, 4, -4
        }
        else
        {
            mtype = 13;
        }
    }
    else // Real
    {
        if ( isSymmetric )
        {
            mtype = 1;
            // also may be: 2, -2
        }
        else
        {
            mtype = 11;
        }
    }

    MKL_INT iparm[64];   /* control parameters */

    for ( int i = 0; i < 64; ++i )
    {
        iparm[i] = 0;
    }

    SCAI_LOG_INFO( logger, "call pardisoinit, mtype = " << mtype )

    pardisoinit( pt, &mtype, iparm );

    for ( int i = 0; i < 64; ++i )
    {
        if ( iparm[i] != 0 )
        {
            SCAI_LOG_DEBUG( logger, "iparm[" << i << "] = " << iparm[i] )
        }
    }

    // iparm has now default values but some changes are required

    if ( ( common::TypeTraits<ValueType>::stype == common::ScalarType::FLOAT ) ||
         ( common::TypeTraits<ValueType>::stype == common::ScalarType::COMPLEX ) )
    {
        iparm[27] = 1;  /* float */
    }
    else if ( ( common::TypeTraits<ValueType>::stype == common::ScalarType::DOUBLE ) ||
              ( common::TypeTraits<ValueType>::stype == common::ScalarType::DOUBLE_COMPLEX ) )
    {
        iparm[27] = 2;  /* double */
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported value type " << common::TypeTraits<ValueType>::stype )
    }

    iparm[34] = 1;  /* PARDISO use C-style indexing for ia and ja arrays */

    MKL_INT phase;
    MKL_INT maxfct = 1;  /* Maximum number of numerical factorizations. */
    MKL_INT mnum   = 1;  /* Which factorization to use. */
    MKL_INT nrhs   = 1;  /* Number of right hand sides */
    MKL_INT msglvl = 0;  /* no statistics */
    MKL_INT error  = 0;  /* Initialize error flag */

    SCAI_LOG_INFO( logger, "Initialization completed ... " )

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */

    phase = 11;

    pardiso( pt, &maxfct, &mnum, &mtype, &phase, const_cast<IndexType*> ( &numRows ),
             const_cast<ValueType*> ( csrValues ), const_cast<IndexType*> ( csrIA ),
             const_cast<IndexType*> ( csrJA ), &iDum, &nrhs, iparm, &msglvl, &vDum, &vDum, &error );

    SCAI_LOG_INFO( logger, "pardiso completed, error = " << error )

    SCAI_PARDISO_ERROR_CHECK ( error, "ERROR during symbolic factorization" )
    SCAI_LOG_INFO( logger, "Reordering completed ... " )
    SCAI_LOG_INFO( logger, "Number of nonzeros in factors = " << iparm[17] );

    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization.                                          */
    /* -------------------------------------------------------------------- */

    phase = 22;

    pardiso( pt, &maxfct, &mnum, &mtype, &phase, const_cast<IndexType*> ( &numRows ),
             const_cast<ValueType*> ( csrValues ), const_cast<IndexType*> ( csrIA ),
             const_cast<IndexType*> ( csrJA ), &iDum, &nrhs, iparm, &msglvl, &vDum, &vDum, &error );

    SCAI_PARDISO_ERROR_CHECK( error, "ERROR during numerical factorization" )
    SCAI_LOG_INFO( logger, "Factorization completed ... " )

    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement.                       */
    /* -------------------------------------------------------------------- */
    phase = 33;
    pardiso( pt, &maxfct, &mnum, &mtype, &phase, const_cast<IndexType*> ( &numRows ),
             const_cast<ValueType*> ( csrValues ), const_cast<IndexType*> ( csrIA ),
             const_cast<IndexType*> ( csrJA ), &iDum, &nrhs, iparm, &msglvl,
             const_cast<ValueType*> ( rhs ), solution, &error );

    SCAI_LOG_INFO( logger, "pardiso 33 completed, error = " << error )
    SCAI_PARDISO_ERROR_CHECK( error, "ERROR during back substitution" )
    SCAI_LOG_INFO( logger, "Solve completed ... " )

    phase = -1; /* Release internal memory. */
    pardiso( pt, &maxfct, &mnum, &mtype, &phase, const_cast<IndexType*> ( &numRows ),
             &vDum, const_cast<IndexType*> ( csrIA ), const_cast<IndexType*> ( csrJA ),
             &iDum, &nrhs, iparm, &msglvl, &vDum, &vDum, &error );

    SCAI_LOG_INFO( logger, "decomposition completed" )
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::ContextType ctx = common::ContextType::Host;
    using kregistry::KernelRegistry;
    SCAI_LOG_INFO( logger, "register CSRUtils MKL-routines for Host at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::decomposition<ValueType> >( decomposition, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

MKLCSRUtils::MKLCSRUtils()
{
    bool useMKL = true;  // by default: use MKL for CSR utils

    common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );

    if ( !useMKL )
    {
        return;
    }

    // replace internal OpenMP kernels

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels( flag );
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

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

MKLCSRUtils MKLCSRUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */

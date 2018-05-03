/**
 * @file OpenMPLAPACK.hpp
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
 * @brief Class with default implementations of LAPACK kernels for host.
 * @author Eric Schricker
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

class COMMON_DLL_IMPORTEXPORT OpenMPLAPACK
{
public:

    /** Implementation of BLASKernelTrait::LAPACK::getrf by LAPACK. */

    template<typename ValueType>
    static void getrf(
        const IndexType m,
        const IndexType n,
        ValueType* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getrf by LAPACK. */

    template<typename ValueType>
    static void getri(
        const IndexType n,
        ValueType* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getinv by LAPACK. */

    template<typename ValueType>
    static void getinv( const IndexType n, ValueType* a, const IndexType lda );

    /** Implementation of BLASKernelTrait::LAPACK::tptrs vi LAPACK. */

    template<typename ValueType>
    static void tptrs(
        const CBLAS_UPLO uplo,
        const common::MatrixOp trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const ValueType* AP,
        ValueType* B,
        const IndexType ldb );

private:

    /** Routine that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPLAPACK();

    /** Destructor for unregistration. */

    ~OpenMPLAPACK();

    /** Static variable for registration at static initialization. */

    static OpenMPLAPACK guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

}; /* OpenMPLAPACK */

} /* end namespace blaskernel */

} /* end namespace scai */

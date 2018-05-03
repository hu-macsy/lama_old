/**
 * @file LAPACK_LAPACK.hpp
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
 * @brief Definition of static class that provides LAPACK kernel implementations via LAPACK.
 * @author Lauretta Schubert
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

class LAPACK_LAPACK
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

    /** Implementation of BLASKernelTrait::LAPACK::getri by LAPACK. */

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
        const common::MatrixOp op,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const ValueType* AP,
        ValueType* B,
        const IndexType ldb );

private:

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

    LAPACK_LAPACK();

    /** Destructor for unregistration. */

    ~LAPACK_LAPACK();

    /** Static variable for registration at static initialization. */

    static LAPACK_LAPACK guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

}; /* LAPACK_LAPACK */

} /* end namespace blaskernel */

} /* end namespace scai */

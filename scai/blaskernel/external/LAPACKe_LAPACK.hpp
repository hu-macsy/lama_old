/**
 * @file LAPACKe_LAPACK.hpp
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
 * @endlicense
 *
 * @brief Implemenation of LAMA LAPACK interface via LAPACKe
 * @author Thomas Brandes
 * @date 10.04.2013
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/OpenMP.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** @brief Implementation of LAPACK functionality at the LAMAInterface by
 *         using LAPACKe (C/C++ extensions for LAPACK).
 *
 *  The use of these routines is more convenient if LAPACKe is available
 *  but has the same functionality as OpenMPLAPACK using Fortran calling conventions.
 */

class LAPACKe_LAPACK
{
public:

    /** Implementation of BLASKernelTrait::LAPACK::getrf by LAPACKe. */

    template<typename ValueType>
    static IndexType getrf(
        const CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        ValueType* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getri by LAPACKe. */

    template<typename ValueType>
    static IndexType getri(
        const CBLAS_ORDER order,
        const IndexType n,
        ValueType* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getinv by LAPACKe. */

    template<typename ValueType>
    static void getinv( const IndexType n, ValueType* a, const IndexType lda );

    /** Implementation of BLASKernelTrait::LAPACK::tptrs by LAPACKe. */

    template<typename ValueType>
    static IndexType tptrs(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const ValueType* AP,
        ValueType* B,
        const IndexType ldb );

private:

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )

    /** Constructor for registration. */

    LAPACKe_LAPACK();

    /** Destructor for unregistration. */

    ~LAPACKe_LAPACK();

    /** Static variable for registration at static initialization. */

    static LAPACKe_LAPACK guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

}; /* LAPACKe_LAPACK */

} /* end namespace blaskernel */

} /* end namespace scai */

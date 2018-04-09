/**
 * @file BLAS_BLAS3.hpp
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
 * @brief BLAS3 utilities on Host Context by wrapping to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>
#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Implementations of blas3 methods for BLASKernelTrait by using BLAS library. */

class COMMON_DLL_IMPORTEXPORT BLAS_BLAS3
{
public:

    /** OpenMP implementation for BLAS3Interface<ValueType>::gemm */

    template<typename ValueType>
    static void gemm(
        const common::MatrixOp opA,
        const common::MatrixOp opB,
        const IndexType M,
        const IndexType N,
        const IndexType K,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* B,
        const IndexType ldb,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

private:

    /** Routine that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    BLAS_BLAS3();

    /** Destructor for unregistration. */

    ~BLAS_BLAS3();

    /** Static variable for registration at static initialization. */

    static BLAS_BLAS3 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */

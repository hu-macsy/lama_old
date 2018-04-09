/**
 * @file OpenMPBLAS2.hpp
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
 * @brief Class with default implementations of BLAS2 routines for host using OpenMP parallelization.
 * @author Eric Schricker
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>
#include <scai/blaskernel/cblas.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** OpenMP Implementations of blas2 methods as specified in BLASKernelTrait.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS2 libraries
 */

class COMMON_DLL_IMPORTEXPORT OpenMPBLAS2
{
public:

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::gemv
     */
    template<typename ValueType>
    static void gemv(
        const common::MatrixOp op,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::geam
     */
    template<typename ValueType>
    static void geam(
        ValueType* C,
        const IndexType ldc,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const common::MatrixOp opA,
        const ValueType beta,
        const ValueType* B,
        const IndexType ldb,
        const common::MatrixOp opB );

private:

    template<typename ValueType>
    static void scale(
        ValueType* C,
        const IndexType ldc,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const common::MatrixOp opA );

    /** structure that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPBLAS2();

    /** Destructor for unregistration. */

    ~OpenMPBLAS2();

    /** Static variable for registration at static initialization. */

    static OpenMPBLAS2 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */

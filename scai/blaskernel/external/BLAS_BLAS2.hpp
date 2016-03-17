/**
 * @file BLAS_BLAS2.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief BLAS_BLAS2.hpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Implementations of blas2 methods of BLASKernelTrait with OpenMP.  */

class COMMON_DLL_IMPORTEXPORT BLAS_BLAS2
{
public:

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::gemv
     */
    template<typename ValueType>
    static void gemv(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE trans,
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

private:

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )

    /** Constructor for registration. */

    BLAS_BLAS2();

    /** Destructor for unregistration. */

    ~BLAS_BLAS2();

    /** Static variable for registration at static initialization. */

    static BLAS_BLAS2 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */

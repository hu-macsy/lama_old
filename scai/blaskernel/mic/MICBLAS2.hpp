/**
 * @file MICBLAS2.hpp
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
 * @brief MICBLAS2.hpp
 * @author Thomas Brandes
 * @date 05.07.2013
 * @since 1.1.0
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

/** Implementations of methods for scai::lama::BLAS2Interface with MIC.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS1 libraries
 */

class COMMON_DLL_IMPORTEXPORT MICBLAS2
{
public:

    /**
     * This function is the MIC implementation of BLASKernelTrait::gemv
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

    SCAI_DECLARE_REGISTRATOR( RegistratorV, template<typename ValueType> )

    /** Helper class for (un) registration of kernel routines at static initialization. */

    class RegisterGuard
    {
    public:
        RegisterGuard();
        ~RegisterGuard();
    };

    static RegisterGuard guard;  // registration of kernels @ static initialization

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */

/**
 * @file OpenMPBLAS3.hpp
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
 * @brief OpenMPBLAS3.hpp
 * @author Eric Schricker
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/blaskernel/cblas.hpp>
#include <scai/logging.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** OpenMP Implementations of blas3 methods as specified in BLASKernelTrait.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS1 libraries
 */

class COMMON_DLL_IMPORTEXPORT OpenMPBLAS3
{
public:

    /** OpenMP implementation for BLASKernelTrait::gemm */

    template<typename ValueType>
    static void gemm(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE TransA,
        const CBLAS_TRANSPOSE TransB,
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

    SCAI_DECLARE_REGISTRATOR( RegistratorV, template<typename ValueType> )

    /** Constructor for registration. */

    OpenMPBLAS3();

    /** Destructor for unregistration. */

    ~OpenMPBLAS3();

    /** Static variable for registration at static initialization. */

    static OpenMPBLAS3 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */

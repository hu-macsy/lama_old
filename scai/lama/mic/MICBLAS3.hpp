/**
 * @file MICBLAS3.hpp
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
 * @brief MICBLAS3.hpp
 * @author Thomas Brandes
 * @date 05.07.2013
 * @since 1.1.0
 */
#ifndef LAMA_MIC_BLAS3_HPP_
#define LAMA_MIC_BLAS3_HPP_

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/LAMATypes.hpp>
#include <scai/lama/SyncToken.hpp>

namespace scai
{

namespace lama
{

/** Implementations of methods for lama::BLAS3Interface with MIC.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS3 libraries
 */

class COMMON_DLL_IMPORTEXPORT MICBLAS3
{
public:

    /** MIC implementation for BLAS3Interface<ValueType>::gemm */

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
        const IndexType ldc,
        SyncToken* syncToken );

    /** Routine that sets functions pointers belonging to BLAS3 in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    static bool initialized;

    static bool registerInterface();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */

#endif // LAMA_MIC_BLAS3_HPP_

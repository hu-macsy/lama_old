/**
 * @file CUDALAPACK.hpp
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
 * @brief CUDALAPACK.hpp
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */
#ifndef SCAI_CUDALAPACK_HPP_
#define SCAI_CUDALAPACK_HPP_

// for dll_import
#include <scai/common/config.hpp>

#include <scai/lama/LAMATypes.hpp>

#include <scai/lama/openmp/BLASHelper.hpp>

namespace tasking
{
    class SyncToken;
}

namespace lama
{


class COMMON_DLL_IMPORTEXPORT CUDALAPACK
{
public:

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    static bool initialized; //!< static initialization used for registration

    static bool registerInterface(); //!< registration

    template<typename ValueType>
    static void laswp(
        const CBLAS_ORDER order,
        const IndexType n,
        ValueType* A,
        const IndexType lda,
        const IndexType k1,
        const IndexType k2,
        const IndexType* ipiv,
        const IndexType incx,
        tasking::SyncToken* syncToken );
};

} /* namespace lama */

#endif // SCAI_CUDALAPACK_HPP_

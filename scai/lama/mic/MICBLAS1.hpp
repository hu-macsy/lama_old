/**
 * @file MICBLAS1.hpp
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
 * @brief MICBLAS1.hpp
 * @author Thomas Brandes
 * @date 05.07.2012
 * @since 1.1.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// logging
#include <scai/logging.hpp>

// others
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace tasking
{
    class SyncToken;
}

namespace lama
{

/** Implementations of methods for scai::lama::BLAS1Interface with MIC.
 *
 *  @todo Add information here about use of native BLAS1 libraries
 */

class COMMON_DLL_IMPORTEXPORT MICBLAS1
{
public:

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::BLAS1::scal
     */
    template<typename ValueType>
    static void scal(
        const IndexType n,
        const ValueType alpha,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::nrm2
     */
    template<typename ValueType>
    static ValueType nrm2( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::asum
     */
    template<typename ValueType>
    static ValueType asum( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::iamax
     */
    template<typename ValueType>
    static IndexType iamax( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::swap
     */
    template<typename ValueType>
    static void swap(
        const IndexType n,
        ValueType* y,
        const IndexType incY,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::copy
     */
    template<typename ValueType>
    static void copy(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::axpy
     */
    template<typename ValueType>
    static void axpy(
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::dot
     */
    template<typename ValueType>
    static ValueType dot(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the MIC implementation of scai::lama::BLAS1Interface::sum
     */
    template<typename ValueType>
    static void sum(
        const IndexType n,
        ValueType alpha,
        const ValueType* x,
        ValueType beta,
        const ValueType* y,
        ValueType* z,
        tasking::SyncToken* syncToken );

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
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

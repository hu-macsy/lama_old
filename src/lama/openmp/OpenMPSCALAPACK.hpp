/**
 * @file OpenMPSCALAPACK.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief OpenMPSCALAPACK.hpp
 * @author lschubert
 * @date 03.07.2012
 * $Id$
 */
#ifndef LAMA_OPENMPSCALAPACK_HPP_
#define LAMA_OPENMPSCALAPACK_HPP_

// others
#include <lama/LAMATypes.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

class OpenMPSCALAPACK
{
public:
    template<typename T>
    static IndexType pgetrf(
        const IndexType m,
        const IndexType n,
        const T* const A,
        const IndexType ia,
        const IndexType ja,
        IndexType* descA,
        IndexType* const ipiv );

    template<typename T>
    static IndexType pgetri(
        const IndexType n,
        const T* const A,
        const IndexType ia,
        const IndexType ja,
        IndexType* const descA,
        IndexType* const ipiv,
        const T* const work,
        IndexType lwork,
        IndexType* const iwork,
        IndexType liwork );

    /** Implementation for SCALAPACKUtils::inverse */

    template<typename T>
    static void inverse( const IndexType n, const IndexType nB, const T* a, const class Communicator& comm );

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

    static bool initialized;

    static bool registerInterface();

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

};

} /* namespace lama */

#endif // LAMA_OPENMPSCALAPACK_HPP_

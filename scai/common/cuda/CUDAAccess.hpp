/**
 * @file CUDAAccess.hpp
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
 * @brief Class to enable/disable access to a CUDA context
 * @author Thomas Brandes
 * @date 04.05.2013
 * @since 1.0.0
 */

#pragma once

#include <cuda.h>

namespace scai
{

namespace common
{

class CUDACtx;

/** This class accesses a CUDA context with the constructor and
 *  releases it with the destructor.
 * 
 *  It guarantees that accesses will also be released in case of exceptions.
 *
 *  Acccesing and releasing a CUDA context is necessary to make it possible
 *  that another thread can also use the same device.
 */

class CUDAAccess
{

public:

    /** The constructor enables the corresponding CUDA context. */

    CUDAAccess( const CUDACtx& dev );

    /** The destructor disables the corresponding CUDA context. */

    ~CUDAAccess();

    /** This method enables an object CUDACtx 
     *
     *  @param[in] ctx context that will be enabled
     *  @returns   pointer to the last enabled context (can be NULL)
     */
    static const CUDACtx* enable( const CUDACtx& ctx );

    /** This method disables the current context and resets the old one.
     *
     *  @param[in] last is pointer to the last context
     */
    static void disable( const CUDACtx* last );

    /** This static method returns the CUDACtx object currently accessed. 
     *
     *  @returns the currently set CUDACtx ( is thread-specific )
     *  @throws Exception if no CUDACtx has been enabled before
     *
     */

    static const CUDACtx& getCurrentCUDACtx();

private:

    CUcontext mCUcontext;

    const CUDACtx* mSaveDevice;  // save the device accessed before
};

} /* end namespace common */

} /* end namespace scai */

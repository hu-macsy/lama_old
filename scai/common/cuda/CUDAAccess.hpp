/**
 * @file CUDAAccess.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Class to enable/disable access to a CUDA context
 * @author Thomas Brandes
 * @date 04.05.2013
 */

#pragma once

#include <scai/common/config.hpp>

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

class COMMON_DLL_IMPORTEXPORT CUDAAccess
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

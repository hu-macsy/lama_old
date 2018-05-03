/**
 * @file CUDACtx.hpp
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
 * @brief Class for a CUDA device
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>


#include <scai/common/NonCopyable.hpp>

#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif ( CUDART_VERSION >= 7050 )
#include <cusolverDn.h>
#include <cusolverSp.h>
#endif

#include <vector>
#include <functional>

namespace scai
{

namespace common
{

/* --------------------------------------------------------------------- */

/** @brief This class defines a context on a CUDA device.
 *
 *  Note: it is possible to define more than one context on a CUDA device,
 *  this is similiar to having multiple processes on one CPU.
 */

class COMMON_DLL_IMPORTEXPORT CUDACtx : private NonCopyable
{

public:

    /** Constructor will check for a CUDA device and create a context for it
     *
     *  @param[in] deviceNr is the number of the device to access
     *  @throws Exception if the device cannot be accessed
     */

    CUDACtx( int deviceNr );

    /** Destructor destroys the CUDA context used for the device.
     *
     *  Be careful: all allocated resources on the device should have been freed before.
     */

    ~CUDACtx();

    /** Getter for the CUcontext to be used for CUDA driver API operations. */

    CUcontext getCUcontext() const
    {
        return mCUcontext;
    }

    /** Getter for the CUdevice to be used for CUDA driver API operations. */

    CUdevice getCUdevice() const
    {
        return mCUdevice;
    }

    /** Getter for the cuBLAS handle, will be created with first use */

    cublasHandle_t getcuBLASHandle() const;

    /** Getter for the cuSparse handle, will be created with first use */

    cusparseHandle_t getcuSparseHandle() const;

#if ( CUDART_VERSION >= 7050 )
    /** Getter for the cuSolverDn handle, will be created with first use */

    cusolverDnHandle_t getcuSolverDnHandle() const;

    /** Getter for the cuSolverSp handle, will be created with first use */

    cusolverSpHandle_t getcuSolverSpHandle() const;
#endif

    /** Getter for the number of the device. */

    int getDeviceNr() const
    {
        return mDeviceNr;
    }

    /** Add a routine to be called with destructor */

    void addShutdown( std::function<void()> routine );

private:

    CUcontext mCUcontext;

    CUdevice mCUdevice;

    int mDeviceNr;

    cublasHandle_t mcuBLASHandle;

    cusparseHandle_t mcuSparseHandle;

#if ( CUDART_VERSION >= 7050 )
    cusolverDnHandle_t mcuSolverDnHandle;
    cusolverSpHandle_t mcuSolverSpHandle;
#endif

    std::vector< std::function<void()> > mShutdownFunctions;
};

}  // namespace common

}  // namespace scai

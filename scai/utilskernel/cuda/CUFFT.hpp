/**
 * @file CUDABLAS1.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief CUDABLAS1.hpp
 * @author lschubert
 * @date 05.07.2012
 */

#pragma once

// local library
#include <scai/utilskernel/FFTKernelTrait.hpp>

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Complex.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

// CUDA
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

namespace scai
{

namespace utilskernel
{

/** Static class that provides CUDA implementaions for the FFT routines as specified in FFTKernelTrait.
 *
 *  The FFT routines are all private and can only be accessed via kernel registry.
 */

class COMMON_DLL_IMPORTEXPORT CUFFT
{

private:

    /** CUDA implementation for FFTKernelTrait::fft */

    template<typename ValueType>
    static void fft(
        common::Complex<ValueType> array[],
        const IndexType nb,
        const IndexType n,
        const IndexType m,
        const int direction );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUFFT();

    /** Destructor for unregistration. */

    ~CUFFT();

    /** Static variable for registration at static initialization. */

    static CUFFT guard;
};

} /* end namespace utilskernel */

} /* end namespace scai */

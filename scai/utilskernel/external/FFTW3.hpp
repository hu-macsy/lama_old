/**
 * @file FFTW3.hpp
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
 * @brief Implementations that uses FFTW3
 * @author Eric Schricker
 * @date 29.09.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT FFTW3
{
public:

    /** Implementation for FFTKernelTrait::fft using FFTW library */

    template<typename ValueType>
    static void fft(
        common::Complex<ValueType> array[],
        const IndexType nb,
        const IndexType n,
        const IndexType m,
        const int direction );

private:

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    FFTW3();

    /** Destructor for unregistration. */

    ~FFTW3();

    /** Static variable for registration at static initialization. */

    static FFTW3 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace utilskernel */

} /* end namespace scai */

/**
 * @file CUDALAPACK.cpp
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
 * @brief Definition of static class with CUDA implementations of LAPACK routines.
 * @author Lauretta Schubert
 * @date 06.07.2012
 */

// hpp
#include <scai/blaskernel/cuda/CUDALAPACK.hpp>

// local library
#include <scai/blaskernel/cuda/CUDABLAS1.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai library
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/cuda/CUDAError.hpp>

namespace scai
{

namespace blaskernel
{

/* ---------------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( CUDALAPACK::logger, "CUDA.LAPACK" )

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDALAPACK::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    // const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_INFO( logger, "register LAPACK routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDALAPACK::CUDALAPACK()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUDALAPACK::~CUDALAPACK()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDALAPACK CUDALAPACK::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */

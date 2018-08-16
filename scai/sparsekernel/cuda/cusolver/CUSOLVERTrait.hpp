/**
 * @file cuda/cusolver/CUSOLVERTrait.hpp
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
 * @brief Definitions for CUSOLVER interface
 * @author Lauretta Schubert
 * @date 19.07.2016
 */

#pragma once

// external
#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif ( CUDART_VERSION >= 7050 )
#include <cusolverDn.h>
#include <cusolverSp.h>

// macros
#define CUSOLVER_DN_NAME( name, prefix ) cusolverDn##prefix##name
#define CUSOLVER_SP_NAME( name, prefix ) cusolverSp##prefix##name

#define CUSOLVER_DN_DEF( name, prefix, retType, definition )          \
    retType CUSOLVER_DN_NAME( name, prefix )( definition );

#define CUSOLVER_SP_DEF( name, prefix, retType, definition )          \
    retType CUSOLVER_SP_NAME( name, prefix )( definition );

#define CUSOLVER_DN_CALL( name, prefix, ... ) \
    SCAI_CUSOLVER_CALL( CUSOLVER_DN_NAME( name, prefix ), __VAR_ARGS__ )

#define CUSOLVER_SP_CALL( name, prefix, ... ) \
    SCAI_CUSOLVER_CALL( CUSOLVER_SP_NAME( name, prefix ), __VAR_ARGS__ )

namespace scai
{

namespace sparsekernel
{

class COMMON_DLL_IMPORTEXPORT CUSOLVERTrait
{
public:
    typedef int BLASIndexType;
    typedef cusolverStatus_t BLASStatus;
    typedef cusolverDnHandle_t SOLVERDnHandle;
    typedef cusolverSpHandle_t SOLVERSpHandle;
    typedef cusparseMatDescr_t BLASMatrix;
};

} /* end namespace sparsekernel */

} /* end namespace scai */

#endif

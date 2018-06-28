/**
 * @file CUBLASTrait.hpp
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
 * @brief Definitions of class to trait calls for methods of the cuBLAS library.
 * @author Eric Stricker
 * @date 21.01.2016
 */

#pragma once

// macros
#define CUBLAS_BLAS_NAME( name, prefix ) cublas##prefix##name

#define CUBLAS_BLAS_DEF( name, prefix, retType, definition )            \
    retType CUBLAS_BLAS_NAME( name, prefix )( definition );

#define CUBLAS_BLAS_CALL( name, prefix, ... )   \
    SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( name, prefix ), __VAR_ARGS__ )

// external
#include <cublas_v2.h>

#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace blaskernel
{

class COMMON_DLL_IMPORTEXPORT CUBLASTrait
{
public:

    typedef int BLASIndexType;
    typedef cublasOperation_t BLASTrans;

    static inline BLASTrans castTrans( const common::MatrixOp op );
};

CUBLASTrait::BLASTrans CUBLASTrait::castTrans( const common::MatrixOp op )
{
    CUBLASTrait::BLASTrans castOp = CUBLAS_OP_N;

    if ( op == common::MatrixOp::NORMAL )
    {
        castOp = CUBLAS_OP_N;
    }
    else if ( op == common::MatrixOp::TRANSPOSE )
    {
        castOp = CUBLAS_OP_T;
    }
    else if ( op == common::MatrixOp::CONJ_TRANSPOSE )
    {
        castOp = CUBLAS_OP_C;
    }

    return castOp;
}

} /* end namespace blaskernel */

} /* end namespace scai */

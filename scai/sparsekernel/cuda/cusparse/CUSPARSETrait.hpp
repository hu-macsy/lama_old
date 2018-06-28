/**
 * @file CUSPARSETrait.hpp
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
 * @brief Definitions for CUBLAS interface
 * @author Eric Stricker
 * @date 21.01.2016
 */

#pragma once

// macros
#define CUSPARSE_BLAS_NAME( name, prefix ) cusparse##prefix##name

#define CUSPARSE_BLAS_DEF( name, prefix, retType, definition )          \
    retType CUSPARSE_BLAS_NAME( name, prefix )( definition );

#define CUSPARSE_BLAS_CALL( name, prefix, ... ) \
    SCAI_CUSPARSE_CALL( CUSPARSE_BLAS_NAME( name, prefix ), __VAR_ARGS__ )

// external
#include <cusparse_v2.h>

namespace scai
{

namespace sparsekernel
{

class COMMON_DLL_IMPORTEXPORT CUSPARSETrait
{
public:
    typedef int BLASIndexType;
    typedef cusparseOperation_t BLASTrans;
    typedef cusparseStatus_t BLASStatus;
    typedef cusparseHandle_t BLASHandle;
    typedef cusparseMatDescr_t BLASMatrix;
    typedef cusparseAction_t BLASOperationType;
    typedef cusparseIndexBase_t BLASIndexBase;
};

} /* end namespace sparsekernel */

} /* end namespace scai */

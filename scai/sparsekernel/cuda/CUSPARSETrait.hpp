/**
 * @file CUBLASTrait.hpp
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
 * @brief Definitions for CUBLAS interface
 * @author  Eric Stricker
 * @date 21.01.2016
 * @since 2.0.0
 */

#pragma once

// macros
#define CUSPARSE_BLAS_NAME( name, prefix ) cusparse##prefix##name

#define CUSPARSE_BLAS_DEF( name, prefix, retType, definition ) 			\
        retType CUSPARSE_BLAS_NAME( name, prefix )( definition );

#define CUSPARSE_BLAS_CALL( name, prefix, ... )	\
		SCAI_CUSPARSE_CALL( CUBLAS_CUSPARSE_NAME( name, prefix ), __VAR_ARGS__ )

// external
#include <cusparse_v2.h>

namespace scai {

namespace sparsekernel {

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

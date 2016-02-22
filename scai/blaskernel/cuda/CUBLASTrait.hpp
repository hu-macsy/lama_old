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
#define CUBLAS_BLAS_NAME( name, prefix ) cublas##prefix##name

#define CUBLAS_BLAS_DEF( name, prefix, retType, definition ) 			\
        retType CUBLAS_BLAS_NAME( name, prefix )( definition );

#define CUBLAS_BLAS_CALL( name, prefix, ... )	\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( name, prefix ), __VAR_ARGS__ )

// external
#include <cublas_v2.h>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT CUBLASTrait
{
public:
	typedef int BLASIndexType;
	typedef cublasOperation_t BLASTrans;
};

} /* end namespace blaskernel */

} /* end namespace scai */

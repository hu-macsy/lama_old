/**
 * @file MICBLASWrapper.hpp
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
 * @endlicense
 *
 * @brief Wrapper for BLAS functions
 * @author Eric Schricker
 * @date 25.01.2016
 */

#pragma once

// local library
#include <scai/blaskernel/cblas.hpp>
#include <scai/blaskernel/mic/MICBLASTrait.hpp>

// scai internal libraries
#include <scai/common/mic/MICCallable.hpp>

// external
#include <mkl_blas.h>

namespace scai {

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MICBLASWrapper;

#define MICBLASWRAPPER_DEF( ValueType, MICBLASValueType, prefix1, prefix2, prefix3, DOT ) 								\
template<>																												\
class COMMON_DLL_IMPORTEXPORT MICBLASWrapper<ValueType>																	\
{																														\
public:																													\
	typedef MICBLASTrait::BLASIndexType BLASIndexType;																	\
	typedef MICBLASTrait::BLASTrans BLASTrans;																			\
																														\
	static MIC_CALLABLE_MEMBER void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,			\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,						\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) 					\
	{																													\
		MIC_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, reinterpret_cast<const MICBLASValueType*>( &alpha ),			\
				reinterpret_cast<const MICBLASValueType*>( A ), &lda, reinterpret_cast<const MICBLASValueType*>( x ), 	\
				&incX, reinterpret_cast<const MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( y ), 	\
				&incY );																								\
	}																													\
																														\
	static MIC_CALLABLE_MEMBER void gemm( const BLASTrans transA, const BLASTrans transB,								\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,										\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,						\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) 						\
	{																													\
		MIC_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, 													\
				reinterpret_cast<const MICBLASValueType*>( &alpha ), reinterpret_cast<const MICBLASValueType*>( A ), 	\
				&lda, reinterpret_cast<const MICBLASValueType*>( B ), &ldb, 											\
				reinterpret_cast<const MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( C ), &ldc );	\
	}																													\
};

MICBLASWRAPPER_DEF( float, float, s, s, is, dot );
MICBLASWRAPPER_DEF( double, double, d, d, id, dot );

#ifdef SCAI_COMPLEX_SUPPORTED
MICBLASWRAPPER_DEF( ComplexFloat, MKL_Complex8, c, sc, ic, dotc );
MICBLASWRAPPER_DEF( ComplexDouble, MKL_Complex16, z, dz, iz, dotc );
#endif

#undef MICBLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */


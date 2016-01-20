/**
 * @file LAPACKeWrapper.hpp
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
 * @brief Wrapper for LAPACKe functions
 * @author Eric Schricker
 * @date 12.11.2015
 * @since 2.0.0
 */

#pragma once


#include <scai/blaskernel/external/LAPACKeDefinitions.hpp>

#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

#include <mkl_lapack.h>
#include <mkl_lapacke.h>

#define FORTRAN_LAPACKE_NAME( name, prefix ) LAPACKE_##prefix##name

namespace scai {

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAPACKeWrapper;

#define LAPACKEWRAPPER_DEF( ValueType, prefix, MKLValueType )							\
template<>																				\
class COMMON_DLL_IMPORTEXPORT LAPACKeWrapper<ValueType>									\
{																						\
public:																					\
	typedef LAPACKeDefinitions::LAPACKIndexType LAPACKIndexType;						\
	typedef LAPACKeDefinitions::LAPACKFlag LAPACKFlag;									\
	typedef LAPACKeDefinitions::LAPACKOrder LAPACKOrder;								\
																						\
	static LAPACKIndexType getrf(const LAPACKOrder matrix_order,						\
			const LAPACKIndexType m, const LAPACKIndexType n,							\
			ValueType* const a, const LAPACKIndexType lda,								\
			LAPACKIndexType* const ipiv)												\
	{																					\
		return FORTRAN_LAPACKE_NAME(getrf, prefix)(matrix_order, m, n, 					\
				reinterpret_cast<MKLValueType*>( a ), lda, ipiv);						\
	}																					\
																						\
	static LAPACKIndexType getri(const LAPACKOrder matrix_order,						\
			const LAPACKIndexType n, ValueType* const A,								\
			const LAPACKIndexType lda,													\
			LAPACKIndexType* const ipiv)												\
	{																					\
		return FORTRAN_LAPACKE_NAME(getri, prefix)(matrix_order, n, 					\
				reinterpret_cast<MKLValueType*>( A ), lda, ipiv);						\
	}																					\
																						\
	static LAPACKIndexType tptrs(const LAPACKOrder matrix_order,						\
			const LAPACKFlag uplo, const LAPACKFlag trans,								\
			const LAPACKFlag diag, const LAPACKIndexType n,								\
			const LAPACKIndexType nrhs, const ValueType* AP,							\
			ValueType* B, const LAPACKIndexType ldb) 									\
	{																					\
		return FORTRAN_LAPACKE_NAME( tptrs, prefix )(matrix_order, uplo, trans, 		\
				diag, n, nrhs, reinterpret_cast<const MKLValueType*>( AP ),				\
				reinterpret_cast<MKLValueType*>( B ), ldb);								\
	}																					\
};

LAPACKEWRAPPER_DEF( float, s, float );
LAPACKEWRAPPER_DEF( double, d, double );
LAPACKEWRAPPER_DEF( ComplexFloat, c, lapack_complex_float );
LAPACKEWRAPPER_DEF( ComplexDouble, z, lapack_complex_double );


} /* end namespace blaskernel */

} /* end namespace scai */


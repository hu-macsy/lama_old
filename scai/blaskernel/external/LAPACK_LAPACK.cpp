/**
 * @file LAPACK_LAPACK.cpp
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
 * @brief LAPACK_LAPACK.cpp
 * @author lschubert
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/external/LAPACK_LAPACK.hpp>

// local library
#include <scai/blaskernel/external/BLAS_BLAS1.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/LAPACKDefinitions.hpp>
#include <scai/blaskernel/external/LAPACKWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

// scai libraries
#include <scai/hmemo/Context.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai {

using common::unique_ptr;

namespace blaskernel {

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER(LAPACK_LAPACK::logger, "LAPACK.LAPACK")

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType LAPACK_LAPACK::getrf(const CBLAS_ORDER order, const IndexType m,
		const IndexType n, ValueType* const A, const IndexType lda,
		IndexType* const ipiv)
{
	SCAI_REGION( "LAPACK.LAPACK.getrf<float>" )

	SCAI_LOG_INFO(logger, "getrf<float> for A of size " << m << " x " << n)

	typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;

	if (common::TypeTraits<IndexType>::stype
			!= common::TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	LAPACKIndexType info = 0;

	if (order == CblasColMajor) {
		info = LAPACKWrapper<ValueType>::getrf(static_cast<LAPACKIndexType>(m),
				static_cast<LAPACKIndexType>(n), A,
				static_cast<LAPACKIndexType>(lda), ipiv);
	} else if (m == n && n == lda) {
		for (IndexType i = 0; i < m; ++i) {
			for (IndexType j = i + 1; j < n; ++j) {
				std::swap(A[i * n + j], A[j * m + i]);
			}
		}

		info = LAPACKWrapper<ValueType>::getrf(static_cast<LAPACKIndexType>(m),
				static_cast<LAPACKIndexType>(n), A,
				static_cast<LAPACKIndexType>(lda), ipiv);

		for (IndexType i = 0; i < m; ++i) {
			for (IndexType j = i + 1; j < n; ++j) {
				std::swap(A[i * n + j], A[j * m + i]);
			}
		}
	} else {
		COMMON_THROWEXCEPTION("row major only supported for square matrices");
	}

	for (IndexType i = 0; i < m; ++i) {
		--ipiv[i]; // Fortran numbering from 1 to n ->  0 to n-1
	}

	if (info < 0) {
		COMMON_THROWEXCEPTION("illegal argument " << ( -info ))
	} else if (info > 0) {
		COMMON_THROWEXCEPTION(
				"value(" << info << "," << info << ")" << " is exactly zero")
	}

	return info;
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::getinv(const IndexType n, ValueType* a,
		const IndexType lda)
{
	SCAI_REGION( "LAPACK.LAPACK.getinv<float>" )

	typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;

	LAPACKIndexType info = 0;

	// unique_ptr, delete by destructor, also done in case of exception

	common::scoped_array<IndexType> ipiv(new IndexType[n]);

	SCAI_LOG_INFO(logger,
			"getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface")

	info = LAPACKWrapper<ValueType>::getrf(static_cast<LAPACKIndexType>(n),
			static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv.get());

	if (info) {
		COMMON_THROWEXCEPTION("LAPACK sgetrf failed, info = " << info)
	}

	common::scoped_array<ValueType> work(new ValueType[n]);

	info = LAPACKWrapper<ValueType>::getri(static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv.get(), work.get(),
			static_cast<LAPACKIndexType>(n));

	if (info) {
		COMMON_THROWEXCEPTION("LAPACK sgetri failed, info = " << info)
	}
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType LAPACK_LAPACK::getri(const CBLAS_ORDER order, const IndexType n,
		ValueType* const a, const IndexType lda, IndexType* const ipiv)
{
	SCAI_REGION( "LAPACK.LAPACK.getri<float>" )

	SCAI_LOG_INFO(logger, "getri<float> for A of size " << n << " x " << n)

	typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;

	if (common::TypeTraits<IndexType>::stype
			!= common::TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	LAPACKIndexType info = 0;

	// translate C indexes into  Fortran Indexes for ipiv

	for (IndexType i = 0; i < n; ++i) {
		++ipiv[i];
	}

	// transpose if not column major order

	if (order != CblasColMajor) {
		SCAI_ASSERT_EQUAL_ERROR(lda, n)

		for (IndexType i = 0; i < n; ++i) {
			// swap row and column

			for (IndexType j = i + 1; j < n; ++j) {
				std::swap(a[i * n + j], a[j * n + i]);
			}
		}
	}

	common::scoped_array<ValueType> work(new ValueType[n]);

	info = LAPACKWrapper<ValueType>::getri(static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv, work.get(),
			static_cast<LAPACKIndexType>(n));

	if (order != CblasColMajor) {
		// transpose back

		for (IndexType i = 0; i < n; ++i) {
			for (IndexType j = i + 1; j < n; ++j) {
				std::swap(a[i * n + j], a[j * n + i]);
			}
		}
	}

	if (info < 0) {
		COMMON_THROWEXCEPTION("illegal argument " << ( -info ))
	} else if (info > 0) {
		COMMON_THROWEXCEPTION(
				"value(" << info << "," << info << ")" << " is exactly zero")
	}

	return info;
}

template<typename ValueType>
IndexType LAPACK_LAPACK::tptrs(const CBLAS_ORDER order, const CBLAS_UPLO uplo,
		const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const IndexType n,
		const IndexType nrhs, const ValueType* AP, ValueType* B,
		const IndexType ldb)
{
	SCAI_REGION( "LAPACK.LAPACK.tptrs<float>" )

	typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;
	typedef LAPACKDefinitions::LAPACKFlag LAPACKFlag;

	LAPACKIndexType info = 0;

	LAPACKFlag UL = LAPACKDefinitions::enum2char(uplo);
	LAPACKFlag TA = LAPACKDefinitions::enum2char(trans);
	LAPACKFlag DI = LAPACKDefinitions::enum2char(diag);

	SCAI_LOG_INFO(logger,
			"tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI);

	if (order == CblasColMajor) {
		LAPACKWrapper<ValueType>::tptrs(UL, TA, DI,
				static_cast<LAPACKIndexType>(n),
				static_cast<LAPACKIndexType>(nrhs), AP, B,
				static_cast<LAPACKIndexType>(ldb));
	} else if (order == CblasRowMajor) {
		// TODO: transpose matrix.
		COMMON_THROWEXCEPTION("row major order not supported for tptrs");
	}

	return info;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::laswp(const CBLAS_ORDER order, const IndexType N,
		ValueType* A, const IndexType LDA, const IndexType K1,
		const IndexType K2, const IndexType* ipiv, const IndexType INCX)
{
	SCAI_REGION( "LAPACK.LAPACK.laswp<float>" )

	typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;

	if (common::TypeTraits<IndexType>::stype
			!= common::TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	if (order == CblasRowMajor) {
		for (IndexType i = K1; i < K2; ++i) {
			if (ipiv[i * INCX] == i) {
				continue;
			}

			BLAS_BLAS1::swap<ValueType>(N, &A[ipiv[i * INCX] * LDA], INCX,
					&A[i * LDA], INCX);
		}
	} else if (order == CblasColMajor) {
		typedef LAPACKDefinitions::LAPACKIndexType LAPACKIndexType;

		LAPACKWrapper<ValueType>::laswp(static_cast<LAPACKIndexType>(N), A,
				static_cast<LAPACKIndexType>(LDA),
				static_cast<LAPACKIndexType>(K1),
				static_cast<LAPACKIndexType>(K2),
				reinterpret_cast<const LAPACKIndexType*>( ipiv ),
				static_cast<LAPACKIndexType>(INCX));
	} else {
		COMMON_THROWEXCEPTION( "cblas_laswp Illegal order setting")
	}
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

void LAPACK_LAPACK::registerKernels(bool deleteFlag)
{
	using kregistry::KernelRegistry;
	using common::context::Host;

	KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE; // priority over OpenMPBLAS

	if (deleteFlag) {
		flag = KernelRegistry::KERNEL_ERASE;
	}

#define LAMA_LAPACK_REGISTER(z, I, _)                                                               \
    KernelRegistry::set<BLASKernelTrait::getrf<ARITHMETIC_HOST_TYPE_##I> >( getrf, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getri<ARITHMETIC_HOST_TYPE_##I> >( getri, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getinv<ARITHMETIC_HOST_TYPE_##I> >( getinv, Host, flag );  \
    KernelRegistry::set<BLASKernelTrait::tptrs<ARITHMETIC_HOST_TYPE_##I> >( tptrs, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::laswp<ARITHMETIC_HOST_TYPE_##I> >( laswp, Host, flag );

	BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_LAPACK_REGISTER, _ )

#undef LAMA_LAPACK_REGISTER
}
	/* --------------------------------------------------------------------------- */
	/*    Static initialiazion at program start                                    */
	/* --------------------------------------------------------------------------- */

LAPACK_LAPACK::LAPACK_LAPACK()
{
	bool deleteFlag = false;
	registerKernels(deleteFlag);
}

LAPACK_LAPACK::~LAPACK_LAPACK()
{
	bool deleteFlag = true;
	registerKernels(deleteFlag);
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

LAPACK_LAPACK LAPACK_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */

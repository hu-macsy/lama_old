/**
 * @file LAPACKe_LAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief LAPACKe_LAPACK.cpp
 * @author Lauretta Schubert
 * @date 02.07.2012
 */

// hpp
#include <scai/blaskernel/external/LAPACKe_LAPACK.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cblas.hpp>
#include <scai/blaskernel/external/LAPACKeWrapper.hpp>
#include <scai/blaskernel/external/LAPACKeTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>

// external
#include <mkl_lapacke.h>

namespace scai 
{

using common::scoped_array;
using common::TypeTraits;

namespace blaskernel 
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER(LAPACKe_LAPACK::logger, "LAPACKe.LAPACK")

/* ------------------------------------------------------------------------- */
/*      getrf                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType LAPACKe_LAPACK::getrf(const CBLAS_ORDER order, const IndexType m,
		const IndexType n, ValueType* const A, const IndexType lda,
		IndexType* const ipiv)
{
	SCAI_LOG_INFO(logger, "getrf<float> for A of size " << m << " x " << n)

	typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

	if (TypeTraits<IndexType>::stype
			!= TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	int info = LAPACKeWrapper<ValueType>::getrf(LAPACKeTrait::enum2order(order),
			static_cast<LAPACKIndexType>(m),
			static_cast<LAPACKIndexType>(n), A,
			static_cast<LAPACKIndexType>(lda), ipiv);

	if (info < 0) {
		COMMON_THROWEXCEPTION("illegal argument " << ( -info ))
	} else if (info > 0) {
		COMMON_THROWEXCEPTION(
				"value(" << info << "," << info << ")" << " is exactly zero")
	}

	return info;
}

/* ------------------------------------------------------------------------- */
/*      getinv                                                               */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::getinv(const IndexType n, ValueType* a,
		const IndexType lda)
{
	typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

	LAPACKIndexType info = 0;

	// scoped_array, will also be freed in case of exception

	scoped_array<LAPACKIndexType> ipiv(
			new LAPACKIndexType[n]);

	SCAI_LOG_INFO(logger,
			"getinv<float> for " << n << " x " << n << " matrix, uses MKL")

	info = LAPACKeWrapper<ValueType>::getrf(LAPACK_COL_MAJOR,
			static_cast<LAPACKIndexType>(n),
			static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv.get());

	// return error if factorization did not work

	if (info) {
		COMMON_THROWEXCEPTION("MKL sgetrf failed, info = " << info)
	}

	info = LAPACKeWrapper<ValueType>::getri(LAPACK_COL_MAJOR,
			static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv.get());

	if (info) {
		COMMON_THROWEXCEPTION("MKL sgetri failed, info = " << info)
	}
}

/* ------------------------------------------------------------------------- */
/*      getri                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int LAPACKe_LAPACK::getri(const CBLAS_ORDER order, const IndexType n,
		ValueType* const a, const IndexType lda, IndexType* const ipiv)
{
	typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

	SCAI_LOG_INFO(logger, "getri<float> for A of size " << n << " x " << n)

	if (TypeTraits<IndexType>::stype
			!= TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	LAPACKeTrait::LAPACKOrder matrix_order = LAPACKeTrait::enum2order(order);

	LAPACKIndexType info = LAPACKeWrapper<ValueType>::getri(matrix_order,
			static_cast<LAPACKIndexType>(n), a,
			static_cast<LAPACKIndexType>(lda), ipiv);

	if (info < 0) {
		COMMON_THROWEXCEPTION("illegal argument " << ( -info ))
	} else if (info > 0) {
		COMMON_THROWEXCEPTION(
				"value(" << info << "," << info << ")" << " is exactly zero")
	}

	return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int LAPACKe_LAPACK::tptrs(const CBLAS_ORDER order, const CBLAS_UPLO uplo,
		const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const IndexType n,
		const IndexType nrhs, const ValueType* AP, ValueType* B,
		const IndexType ldb)
{
	typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

	LAPACKeTrait::LAPACKFlag UL = LAPACKeTrait::enum2char(uplo);
	LAPACKeTrait::LAPACKFlag TA = LAPACKeTrait::enum2char(trans);
	LAPACKeTrait::LAPACKFlag DI = LAPACKeTrait::enum2char(diag);

	LAPACKeTrait::LAPACKOrder matrix_order = LAPACKeTrait::enum2order(order);

	if (TypeTraits<IndexType>::stype
			!= TypeTraits<LAPACKIndexType>::stype) {
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION("indextype mismatch");
	}

	SCAI_LOG_INFO(logger,
			"tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI);

	SCAI_ASSERT_ERROR(ldb >= std::max(1, n), "ldb = " << ldb << " out of range");

	int info = LAPACKeWrapper<ValueType>::tptrs(matrix_order, UL, TA, DI,
			static_cast<LAPACKIndexType>(n),
			static_cast<LAPACKIndexType>(nrhs), AP, B,
			static_cast<LAPACKIndexType>(ldb));

	return info;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
	using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register lapack wrapper routines for Host at kernel registry" )

    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( LAPACKe_LAPACK::getrf, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( LAPACKe_LAPACK::getri, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( LAPACKe_LAPACK::getinv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( LAPACKe_LAPACK::tptrs, ctx, flag );
}

LAPACKe_LAPACK::LAPACKe_LAPACK() {
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_REPLACE );
}

LAPACKe_LAPACK::~LAPACKe_LAPACK() {
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

LAPACKe_LAPACK LAPACKe_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */

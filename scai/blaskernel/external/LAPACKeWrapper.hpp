/**
 * @file LAPACKeWrapper.hpp
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
 * @brief Definition of wrapper class for LAPACKe functions
 * @author Eric Schricker
 * @date 12.11.2015
 */

#pragma once

// local library
#include <scai/blaskernel/external/LAPACKeTrait.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

// external
#include <mkl_lapacke.h>

namespace scai
{

namespace blaskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAPACKeWrapper;

#define LAPACKEWRAPPER_DEF( ValueType, MKLValueType, prefix )                           \
    template<>                                                                              \
    class COMMON_DLL_IMPORTEXPORT LAPACKeWrapper<ValueType>                                 \
    {                                                                                       \
    public:                                                                                 \
        typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;                              \
        typedef LAPACKeTrait::LAPACKFlag LAPACKFlag;                                        \
        typedef LAPACKeTrait::LAPACKOrder LAPACKOrder;                                      \
        \
        static LAPACKIndexType getrf(const LAPACKOrder matrix_order,                        \
                                     const LAPACKIndexType m, const LAPACKIndexType n,                           \
                                     ValueType* const a, const LAPACKIndexType lda,                              \
                                     LAPACKIndexType* const ipiv)                                                \
        {                                                                                   \
            return FORTRAN_LAPACKE_NAME(getrf, prefix)(matrix_order, m, n,                  \
                    reinterpret_cast<MKLValueType*>( a ), lda, ipiv);                       \
        }                                                                                   \
        \
        static LAPACKIndexType getri(const LAPACKOrder matrix_order,                        \
                                     const LAPACKIndexType n, ValueType* const A,                                \
                                     const LAPACKIndexType lda,                                                  \
                                     LAPACKIndexType* const ipiv)                                                \
        {                                                                                   \
            return FORTRAN_LAPACKE_NAME(getri, prefix)(matrix_order, n,                     \
                    reinterpret_cast<MKLValueType*>( A ), lda, ipiv);                       \
        }                                                                                   \
        \
        static LAPACKIndexType tptrs(const LAPACKOrder matrix_order,                        \
                                     const LAPACKFlag uplo, const LAPACKFlag trans,                              \
                                     const LAPACKFlag diag, const LAPACKIndexType n,                             \
                                     const LAPACKIndexType nrhs, const ValueType* AP,                            \
                                     ValueType* B, const LAPACKIndexType ldb)                                    \
        {                                                                                   \
            return FORTRAN_LAPACKE_NAME( tptrs, prefix )(matrix_order, uplo, trans,         \
                    diag, n, nrhs, reinterpret_cast<const MKLValueType*>( AP ),             \
                    reinterpret_cast<MKLValueType*>( B ), ldb);                             \
        }                                                                                   \
    };

LAPACKEWRAPPER_DEF( float, float, s );
LAPACKEWRAPPER_DEF( double, double, d );

#ifdef SCAI_COMPLEX_SUPPORTED
LAPACKEWRAPPER_DEF( ComplexFloat, lapack_complex_float, c );
LAPACKEWRAPPER_DEF( ComplexDouble, lapack_complex_double, z );
#endif

#undef LAPACKEWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */


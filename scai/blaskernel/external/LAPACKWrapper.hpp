/**
 * @file LAPACKWrapper.hpp
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
 * @brief Wrapper for LAPACK functions
 * @author Eric Schricker
 * @date 14.01.2016
 */

#pragma once

// local library
#include <scai/blaskernel/external/LAPACKTrait.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace blaskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAPACKWrapper;

#define LAPACKWRAPPER_DEF( ValueType, prefix )                                                          \
    template<>                                                                                              \
    class COMMON_DLL_IMPORTEXPORT LAPACKWrapper<ValueType>                                                  \
    {                                                                                                       \
    public:                                                                                                 \
        typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;                                               \
        typedef LAPACKTrait::LAPACKFlag LAPACKFlag;                                                         \
        \
        static LAPACKIndexType getrf(                                                                       \
                const LAPACKIndexType m,                                                                    \
                const LAPACKIndexType n, ValueType* a,                                                      \
                const LAPACKIndexType lda,                                                                  \
                LAPACKIndexType* ipivot)                                                                    \
        {                                                                                                   \
            LAPACKIndexType info;                                                                           \
            FORTRAN_LAPACK_NAME( getrf, prefix )(&m, &n, a, &lda, ipivot, &info);                           \
            return info;                                                                                    \
        }                                                                                                   \
        \
        static LAPACKIndexType getri(                                                                       \
                const LAPACKIndexType n, ValueType* a,                                                      \
                const LAPACKIndexType lda,                                                                  \
                LAPACKIndexType* ipivot, ValueType* work,                                                   \
                const LAPACKIndexType ldwork)                                                               \
        {                                                                                                   \
            LAPACKIndexType info;                                                                           \
            FORTRAN_LAPACK_NAME( getri, prefix )(&n, a, &lda, ipivot, work, &ldwork, &info);                \
            return info;                                                                                    \
        }                                                                                                   \
        \
        static LAPACKIndexType tptrs(LAPACKFlag uplo,                                                       \
                                     LAPACKFlag transa, LAPACKFlag diag,                                    \
                                     const LAPACKIndexType n,                                               \
                                     const LAPACKIndexType nrhs, const ValueType* ap,                       \
                                     ValueType* b, const LAPACKIndexType ldb)                               \
        {                                                                                                   \
            LAPACKIndexType info;                                                                           \
            FORTRAN_LAPACK_NAME( tptrs, prefix )(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );    \
            return info;                                                                                    \
        }                                                                                                   \
    };


LAPACKWRAPPER_DEF( float, s )
LAPACKWRAPPER_DEF( double, d )

#ifdef SCAI_COMPLEX_SUPPORTED
LAPACKWRAPPER_DEF( ComplexFloat, c )
LAPACKWRAPPER_DEF( ComplexDouble, z )
#endif

#undef LAPACKWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */


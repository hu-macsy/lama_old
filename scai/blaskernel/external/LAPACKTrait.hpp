/**
 * @file LAPACKTrait.hpp
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
 * @brief Trait of LAPACK functions
 * @author Eric Schricker
 * @date 14.01.2016
 */

#pragma once

// internal scai libraries
#include <scai/common/config.hpp>
#include <scai/common/MatrixOp.hpp>

//macros
#define FORTRAN_LAPACK_NAME( name, prefix ) prefix##name##_

#define FORTRAN_LAPACK_DEF( name, prefix, retType, definition )             \
    retType FORTRAN_LAPACK_NAME( name, prefix )( definition );

namespace scai
{

namespace blaskernel
{

class COMMON_DLL_IMPORTEXPORT LAPACKTrait
{
public:
#ifdef F77_INT
    typedef F77_INT LAPACKIndexType;
#else
    typedef int LAPACKIndexType;
#endif

#ifdef F77_CHAR
    typedef F77_CHAR LAPACKFlag;
#else
    typedef char LAPACKFlag;
#endif

    static inline LAPACKFlag enum2char( const CBLAS_UPLO uplo )
    {
        switch ( uplo )
        {
            case CblasUpper:
                return 'U';

            case CblasLower:
                return 'L';

            default:
                COMMON_THROWEXCEPTION( "Illegal uplo: " << uplo );
        }
    }

    static inline LAPACKFlag enum2char( const common::MatrixOp op )
    {
        switch ( op )
        {
            case common::MatrixOp::NORMAL:
                return 'N';

            case common::MatrixOp::TRANSPOSE:
                return 'T';

            case common::MatrixOp::CONJ_TRANSPOSE:
                return 'C';

            default:
                COMMON_THROWEXCEPTION( "Illegal matrix op: " << op );
        }
    }

    static inline LAPACKFlag enum2char( const CBLAS_DIAG diag )
    {
        switch ( diag )
        {
            case CblasNonUnit:
                return 'N';

            case CblasUnit:
                return 'U';

            default:
                COMMON_THROWEXCEPTION( "Illegal diag: " << diag );
        }
    }
};

} /* end namespace blaskernel */

} /* end namespace scai */

extern "C"
{
#define CALL_DEF_GETRI( ValueType ) const scai::blaskernel::LAPACKTrait::LAPACKIndexType* n, ValueType* a, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* lda, scai::blaskernel::LAPACKTrait::LAPACKIndexType* ipivot, ValueType* work, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* ldwork, scai::blaskernel::LAPACKTrait::LAPACKIndexType* info
#define CALL_DEF_GETRF( ValueType ) const scai::blaskernel::LAPACKTrait::LAPACKIndexType* m,    const scai::blaskernel::LAPACKTrait::LAPACKIndexType* n, ValueType* a,const scai::blaskernel::LAPACKTrait::LAPACKIndexType* lda, scai::blaskernel::LAPACKTrait::LAPACKIndexType* ipivot, scai::blaskernel::LAPACKTrait::LAPACKIndexType* info
#define CALL_DEF_TPTRS( ValueType ) const scai::blaskernel::LAPACKTrait::LAPACKFlag* uplo, const scai::blaskernel::LAPACKTrait::LAPACKFlag* transa, const scai::blaskernel::LAPACKTrait::LAPACKFlag* diag, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* n, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* nrhs, const ValueType* ap, ValueType* b, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* ldb, scai::blaskernel::LAPACKTrait::LAPACKIndexType* info
#define CALL_DEF_LASWP( ValueType ) const scai::blaskernel::LAPACKTrait::LAPACKIndexType* n, ValueType* a, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* lda, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* k1, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* k2, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* ipiv, const scai::blaskernel::LAPACKTrait::LAPACKIndexType* incx

// float
    FORTRAN_LAPACK_DEF( getri, s, void, CALL_DEF_GETRI( float ) )
    FORTRAN_LAPACK_DEF( getrf, s, void, CALL_DEF_GETRF( float ) )
    FORTRAN_LAPACK_DEF( tptrs, s, void, CALL_DEF_TPTRS( float ) )
// double
    FORTRAN_LAPACK_DEF( getrf, d, void, CALL_DEF_GETRF( double ) )
    FORTRAN_LAPACK_DEF( getri, d, void, CALL_DEF_GETRI( double ) )
    FORTRAN_LAPACK_DEF( tptrs, d, void, CALL_DEF_TPTRS( double ) )

#ifdef SCAI_COMPLEX_SUPPORTED
// ComplexFloat
    FORTRAN_LAPACK_DEF( getrf, c, void, CALL_DEF_GETRF( scai::ComplexFloat ) )
    FORTRAN_LAPACK_DEF( getri, c, void, CALL_DEF_GETRI( scai::ComplexFloat ) )
    FORTRAN_LAPACK_DEF( tptrs, c, void, CALL_DEF_TPTRS( scai::ComplexFloat ) )
// ComplexDouble
    FORTRAN_LAPACK_DEF( getri, z, void, CALL_DEF_GETRI( scai::ComplexDouble ) )
    FORTRAN_LAPACK_DEF( getrf, z, void, CALL_DEF_GETRF( scai::ComplexDouble ) )
    FORTRAN_LAPACK_DEF( tptrs, z, void, CALL_DEF_TPTRS( scai::ComplexDouble ) )
#endif

#undef CALL_DEF_GETRI
#undef CALL_DEF_GETRF
#undef CALL_DEF_TPTRS
#undef CALL_DEF_LASWP

#undef FORTRAN_LAPACK_DEF

} /* end extern "C" */

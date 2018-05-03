/**
 * @file LAPACKeTrait.hpp
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
 * @brief Typedefs and macros for wrapping routines of the LAPACKe library.
 * @author Eric Schricker
 * @date 14.01.2016
 */

#pragma once

// internal scai libraries
#include <scai/common/config.hpp>

// external
#include <mkl_lapacke.h>

// macros
#define FORTRAN_LAPACKE_NAME( name, prefix ) LAPACKE_##prefix##name

namespace scai
{

namespace blaskernel
{

class COMMON_DLL_IMPORTEXPORT LAPACKeTrait
{
public:
    typedef lapack_int LAPACKIndexType;
    typedef char LAPACKFlag;
    typedef int LAPACKOrder;

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

    static inline LAPACKOrder enum2order( const CBLAS_ORDER order )
    {
        switch ( order )
        {
            case CblasColMajor:
                return LAPACK_COL_MAJOR;

            case CblasRowMajor:
                return LAPACK_ROW_MAJOR;

            default:
                COMMON_THROWEXCEPTION( "illegal matrix order " << order )
        }
    }
};

} /* end namespace blaskernel */

} /* end namespace scai */


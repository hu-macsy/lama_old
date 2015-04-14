/**
 * @file BLASHelper.cpp
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
 * @brief BLASHelper.cpp
 * @author Lauretta Schubert
 * @date 10.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/BLASHelper.hpp>

// macros
#include <lama/macros/unused.hpp>

#include <sstream>

namespace lama
{

/* ------------------------------------------------------------------------- */

char BLASHelper::lapack_uplo( const CBLAS_UPLO uplo )
{
    char UL = 'U';

    if( uplo == CblasUpper )
    {
        UL = 'U';
    }
    else if( uplo == CblasLower )
    {
        UL = 'L';
    }
    else
    {
        LAMA_THROWEXCEPTION( "Illegal uplo: " << uplo );
    }

    return UL;
}

/* ------------------------------------------------------------------------- */

char BLASHelper::lapack_transpose( const CBLAS_TRANSPOSE trans )
{
    char TA = 'N';

    if( trans == CblasNoTrans )
    {
        TA = 'N';
    }
    else if( trans == CblasTrans )
    {
        TA = 'T';
    }
    else if( trans == CblasConjTrans )
    {
        TA = 'C';
    }
    else
    {
        LAMA_THROWEXCEPTION( "Illegal trans: " << trans );
    }

    return TA;
}

/* ------------------------------------------------------------------------- */

char BLASHelper::lapack_diag( const CBLAS_DIAG diag )
{
    char DI = 'N';

    if( diag == CblasNonUnit )
    {
        DI = 'N';
    }
    else if( diag == CblasUnit )
    {
        DI = 'U';
    }
    else
    {
        LAMA_THROWEXCEPTION( "Illegal diag: " << diag );
    }

    return DI;
}

/* ------------------------------------------------------------------------- */

void BLASHelper::XERBLA_cpu(
    int UNUSED(RowMajorStrg),
    int UNUSED(info),
    const char *UNUSED(rout),
    const char *UNUSED(form),
    ... )
{
// @todo This routine should proably call blas error function xerbla ?!
}

} /* namespace lama */

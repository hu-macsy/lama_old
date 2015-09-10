/**
 * @file BLASHelper.hpp
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
 * @brief BLASHelper.hpp
 * @author Lauretta Schubert
 * @date 04.07.2012
 * @since 1.0.0
 */

#pragma once

#include <scai/lama/macros/inline.hpp>
#include <scai/lama/LAMATypes.hpp>
#include <scai/common/exception/Exception.hpp>

#include <scai/lama/cblas.hpp>

namespace scai
{

namespace lama
{

class BLASHelper
{
public:

    static char lapack_transpose( const CBLAS_TRANSPOSE trans );
    static char lapack_uplo( const CBLAS_UPLO uplo );
    static char lapack_diag( const CBLAS_DIAG diag );

    void cblas_xerbla( int p, const char *rout, const char *form, ... );
    //    void setComputeUnit(enum COMPUTE_UNIT unit);
    //    enum COMPUTE_UNIT getComputeUnit();

    static void XERBLA_cpu( int RowMajorStrg, int info, const char *rout, const char *form, ... );
};

/*
 * Converts a floatPointerPointer to a voidPointerPointer without emmiting a Compiler warning.
 * Results in the same operation as casting to a void Pointer( tested on gcc 4.6.0(20110429) with -O3).
 * There is a little Overhead in the unoptimized version. But that shouldn't matter in a debug build.
 */
LAMA_STATIC_INLINE_FUNCTION_PREFIX void** lama_sToVoidPtr( float** floatPointer )
{
    union
    {
        float** fpp;
        void** vpp;
    } convert;
    convert.fpp = floatPointer;
    return convert.vpp;
}

/*
 * Converts a doublePointerPointer to a voidPointerPointer without emmiting a Compiler warning.
 * Results in the same operation as casting to a void Pointer( tested on gcc 4.6.0(20110429) with -O3).
 * There is a little Overhead in the unoptimized version. But that shouldn't matter in a debug build.
 */LAMA_STATIC_INLINE_FUNCTION_PREFIX void** lama_dToVoidPtr( double** doublePointer )
{
    union
    {
        double** dpp;
        void** vpp;
    } convert;
    convert.dpp = doublePointer;
    return convert.vpp;
}

} /* end namespace lama */

} /* end namespace scai */

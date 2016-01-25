/**
 * @file cublas_cast.hpp
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
 * @brief Help routines for interface from CUDA Kernels to cuBLAS
 * @author Thomas Brandes
 * @date 05.06.2014
 * @since 1.1.0
 */

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

// CUDA
#include <cublas_v2.h>

namespace scai
{

namespace blaskernel
{

#ifdef SCAI_COMPLEX_SUPPORTED

/* ---------------------------------------------------------------------------------------*/
/*    cublasCast                                                                          */
/* ---------------------------------------------------------------------------------------*/

/*  cublasCast converts pointers to LAMA complex numbers to
 *  cuBlas pointers for complex numbers. This is safe as both
 *  are internally represented in the same way.
 *
 *  Note: these routines are always called on the __host__ and never on __device__
 */

/**
 * @brief convert pointer to ComplexFloat to pointer cuComplex
 */
static inline cuComplex* cublasCast( ComplexFloat* x )
{
    return reinterpret_cast<cuComplex*>( x );
}

/**
 * @brief convert pointer to ComplexDouble to pointer cuDoubleComplex
 */
static inline cuDoubleComplex* cublasCast( ComplexDouble* x )
{
    return reinterpret_cast<cuDoubleComplex*>( x );
}

/**
 * @brief convert const pointer to ComplexFloat to const pointer cuComplex
 */
static inline const cuComplex* cublasCast( const ComplexFloat* x )
{
    return reinterpret_cast<const cuComplex*>( x );
}

/**
 * @brief convert const pointer to ComplexDouble to const pointer cuDoubleComplex
 */
static inline const cuDoubleComplex* cublasCast( const ComplexDouble* x )
{
    return reinterpret_cast<const cuDoubleComplex*>( x );
}

/**
 * @brief convert value ComplexFloat to value cuComplex
 */
static inline cuComplex cublasCast( ComplexFloat x )
{
    return *cublasCast( &x );
}

/**
 * @brief convert value ComplexDouble to value cuDoubleComplex
 */
static inline cuDoubleComplex cublasCast( ComplexDouble x )
{
    return *cublasCast( &x );
}

#endif

} /* end namespace blaskernel */

} /* end namespace scai */

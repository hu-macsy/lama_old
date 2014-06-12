/**
 * @file cblas.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief C Inteface to BLAS routines
 * @author Thomas Brandes
 * @date 05.06.2014
 * @since 1.1.0
 */

#ifndef LAMA_CBLAS_HPP_
#define LAMA_CBLAS_HPP_

/* Note: LAMA uses these type definitions:

#define CBLAS_INDEX size_t  

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

*/

/** Solution 1: use cblas interface of MKL
 *  
 *  - recommended when using MKL
 *  - mandatory when also using other MKL routines (e.g. for CSR matrices)
 */

#include <mkl_cblas.h>

//  #include <lama/cblas/include/cblas.h>

#endif //  LAMA_CBLAS_HPP_

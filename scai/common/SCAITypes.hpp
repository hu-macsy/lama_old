/**
 * @file SCAITypes.hpp
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
 * @brief SCAITypes.hpp
 * @author Jiri Kraus
 * @date 23.02.2011
 */
#pragma once

// local library

// no support of Complex if this file is not included

#include <scai/common/Complex.hpp>

// std
#include <cstring>
#include <limits>
#include <stdint.h>

/** Common namespace for all projects of Fraunhofer SCAI. */

namespace scai
{
}

/** LAMA uses for all its classes and routines an own namespace.
 *
 *  Applications using LAMA must either put \c scai::lama:: before the used
 *  classes and routines or utilize the using directive of C++.
 *
 *  \code
 *     scai::lama::DenseVector<float> V;
 *     using namespace scai::lama;
 *     DenseVector<float> V;
 *  \endcode
 */

/** Data type that is used for indexing in LAMA.
 *
 *  Currently, it is still necessary that it is a signed data type.
 *  int is the good choice, might be long int or long long int for
 *  future versions that deal with very large matrices even on on processor.
 */
typedef int IndexType;

/** Data type long double to have it as one word. Otherwise certain macros
 *  might fail to work correctly.
 */
typedef long double LongDouble;

/** Definition for a constant value that indicates a non-available index.
 */

#ifdef __INTEL_OFFLOAD
__declspec( target(mic) )
#endif
extern const IndexType nIndex;

/** Data type that is used for numbering of partitions.
 *
 */
typedef IndexType PartitionId;

/** Definition for a constant value that indicates a non-available partition.
 */
#ifdef __INTEL_OFFLOAD
__declspec( target(mic) )
#endif
extern const PartitionId nPartition;

// List here all arithmetic types for which matrices, storages might be created

#define ARITHMETIC_HOST_TYPE_0 float
#define ARITHMETIC_HOST_TYPE_1 double
#define ARITHMETIC_HOST_TYPE_2 ComplexFloat
#define ARITHMETIC_HOST_TYPE_3 ComplexDouble
#define ARITHMETIC_HOST_TYPE_4 long double
#define ARITHMETIC_HOST_TYPE_5 ComplexLongDouble

/** Number of supported types used in REPEAT macros */
#ifdef SCAI_COMPLEX_SUPPORTED

#define ARITHMETIC_HOST_TYPE_CNT 6

#else

#define ARITHMETIC_HOST_TYPE_CNT 2

#endif

/** Number of supported types in external libraries like BLAS, LAPACK */
#ifdef SCAI_COMPLEX_SUPPORTED

#define ARITHMETIC_HOST_EXT_TYPE_CNT 4

#else

#define ARITHMETIC_HOST_EXT_TYPE_CNT 2

#endif

// List CUDA values on its own

#define ARITHMETIC_CUDA_TYPE_0 float
#define ARITHMETIC_CUDA_TYPE_1 double
#define ARITHMETIC_CUDA_TYPE_2 ComplexFloat
#define ARITHMETIC_CUDA_TYPE_3 ComplexDouble

/** Number of suported types by CUDA devices */
#ifdef SCAI_COMPLEX_SUPPORTED

#define ARITHMETIC_CUDA_TYPE_CNT 4

#else

#define ARITHMETIC_CUDA_TYPE_CNT 2

#endif


// List MIC value types on its own

#define ARITHMETIC_MIC_TYPE_0 float
#define ARITHMETIC_MIC_TYPE_1 double
#define ARITHMETIC_MIC_TYPE_2 ComplexFloat
#define ARITHMETIC_MIC_TYPE_3 ComplexDouble

/** Number of supported types by MIC devices */
#ifdef SCAI_COMPLEX_SUPPORTED

#define ARITHMETIC_MIC_TYPE_CNT 4

#else

#define ARITHMETIC_MIC_TYPE_CNT 2

#endif

// For convenience we define ARRAY_TYPE, must be ARITHMETIC_HOST_TYPE_CNT + 1


#define ARRAY_TYPE0    IndexType
#define ARRAY_TYPE1    float
#define ARRAY_TYPE2    double
#define ARRAY_TYPE3    ComplexFloat
#define ARRAY_TYPE4    ComplexDouble
#define ARRAY_TYPE5    long double
#define ARRAY_TYPE6    ComplexLongDouble

#ifdef SCAI_COMPLEX_SUPPORTED

#define ARRAY_TYPE_CNT 7

#else

#define ARRAY_TYPE_CNT 3

#endif

#ifdef SCAI_COMPLEX_SUPPORTED

typedef ComplexLongDouble ScalarRepType;

#else

typedef double ScalarRepType;

#endif

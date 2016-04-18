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
 * @brief Define all arithmetic types for which SCAI template classes will be instantiated
 * @author Jiri Kraus
 * @date 23.02.2011
 */
#pragma once

// local library

// no support of Complex if this file is not included

#include <scai/common/Complex.hpp>
#include <scai/common/mic/MICCallable.hpp>
#include <scai/common/mepr/TypeList.hpp>

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

MIC_CALLABLE_MEMBER extern const IndexType nIndex;

/** Data type that is used for numbering of partitions.
 *
 */
typedef IndexType PartitionId;

/** Definition for a constant value that indicates a non-available partition.
 */
MIC_CALLABLE_MEMBER extern const PartitionId nPartition;

// List here all arithmetic types for which matrices, storages might be created

#define SCAI_ARITHMETIC_HOST_TYPE_0 float
#define SCAI_ARITHMETIC_HOST_TYPE_1 double
#define SCAI_ARITHMETIC_HOST_TYPE_2 ComplexFloat
#define SCAI_ARITHMETIC_HOST_TYPE_3 ComplexDouble
#define SCAI_ARITHMETIC_HOST_TYPE_4 long double
#define SCAI_ARITHMETIC_HOST_TYPE_5 ComplexLongDouble

/*
 * Kernel
 */

#define SCAI_ARITHMETIC_HOST_CNT 6
#define SCAI_ARITHMETIC_HOST     float, double, ComplexFloat, ComplexDouble, long double, ComplexLongDouble

#define SCAI_ARITHMETIC_EXT_HOST_CNT 4
#define SCAI_ARITHMETIC_EXT_HOST float, double, ComplexFloat, ComplexDouble

#define SCAI_ARITHMETIC_CUDA_CNT 4
#define SCAI_ARITHMETIC_CUDA     float, double, ComplexFloat, ComplexDouble

#define SCAI_ARITHMETIC_MIC_CNT 2
#define SCAI_ARITHMETIC_MIC      float, double

#define SCAI_ARITHMETIC_HOST_LIST TYPELIST(     SCAI_ARITHMETIC_HOST_CNT,     SCAI_ARITHMETIC_HOST )
#define SCAI_ARITHMETIC_EXT_HOST_LIST TYPELIST( SCAI_ARITHMETIC_EXT_HOST_CNT, SCAI_ARITHMETIC_EXT_HOST )
#define SCAI_ARITHMETIC_CUDA_LIST TYPELIST(     SCAI_ARITHMETIC_CUDA_CNT,     SCAI_ARITHMETIC_CUDA )
#define SCAI_ARITHMETIC_MIC_LIST TYPELIST(      SCAI_ARITHMETIC_MIC_CNT,      SCAI_ARITHMETIC_MIC )

/*
 * Array
 */

#define SCAI_ARITHMETIC_ARRAY_HOST_CNT 7
#define SCAI_ARITHMETIC_ARRAY_HOST     IndexType, SCAI_ARITHMETIC_HOST

#define SCAI_ARITHMETIC_ARRAY_EXT_HOST_CNT 5
#define SCAI_ARITHMETIC_ARRAY_EXT_HOST IndexType, SCAI_ARITHMETIC_EXT_HOST

#define SCAI_ARITHMETIC_ARRAY_CUDA_CNT 5
#define SCAI_ARITHMETIC_ARRAY_CUDA     IndexType, SCAI_ARITHMETIC_CUDA

#define SCAI_ARITHMETIC_ARRAY_MIC_CNT 3
#define SCAI_ARITHMETIC_ARRAY_MIC      IndexType, SCAI_ARITHMETIC_MIC

#define SCAI_ARITHMETIC_ARRAY_HOST_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_HOST_CNT, SCAI_ARITHMETIC_ARRAY_HOST )
#define SCAI_ARITHMETIC_ARRAY_CUDA_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_CUDA_CNT, SCAI_ARITHMETIC_ARRAY_CUDA )
#define SCAI_ARITHMETIC_ARRAY_MIC_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_MIC_CNT, SCAI_ARITHMETIC_ARRAY_MIC )



/** Number of supported types used in REPEAT macros */
#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_ARITHMETIC_HOST_TYPE_CNT 6

#else

#define SCAI_ARITHMETIC_HOST_TYPE_CNT 2

#endif

/** Number of supported types in external libraries like BLAS, LAPACK */
#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_ARITHMETIC_HOST_EXT_TYPE_CNT 4

#else

#define SCAI_ARITHMETIC_HOST_EXT_TYPE_CNT 2

#endif

// List CUDA values on its own

#define SCAI_ARITHMETIC_CUDA_TYPE_0 float
#define SCAI_ARITHMETIC_CUDA_TYPE_1 double
#define SCAI_ARITHMETIC_CUDA_TYPE_2 ComplexFloat
#define SCAI_ARITHMETIC_CUDA_TYPE_3 ComplexDouble

/** Number of suported types by CUDA devices */
#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_ARITHMETIC_CUDA_TYPE_CNT 4

#else

#define SCAI_ARITHMETIC_CUDA_TYPE_CNT 2

#endif


// List MIC value types on its own

#define SCAI_ARITHMETIC_MIC_TYPE_0 float
#define SCAI_ARITHMETIC_MIC_TYPE_1 double
#define SCAI_ARITHMETIC_MIC_TYPE_2 ComplexFloat
#define SCAI_ARITHMETIC_MIC_TYPE_3 ComplexDouble

/** Number of supported types by MIC devices */
#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_ARITHMETIC_MIC_TYPE_CNT 4

#else

#define SCAI_ARITHMETIC_MIC_TYPE_CNT 2

#endif

// For convenience we define SCAI_ARRAY_TYPE, must be SCAI_ARITHMETIC_HOST_TYPE_CNT + 1


#define SCAI_ARRAY_TYPE0    IndexType
#define SCAI_ARRAY_TYPE1    float
#define SCAI_ARRAY_TYPE2    double
#define SCAI_ARRAY_TYPE3    ComplexFloat
#define SCAI_ARRAY_TYPE4    ComplexDouble
#define SCAI_ARRAY_TYPE5    long double
#define SCAI_ARRAY_TYPE6    ComplexLongDouble

#ifdef SCAI_COMPLEX_SUPPORTED

#define SCAI_ARRAY_TYPE_CNT 7

#else

#define SCAI_ARRAY_TYPE_CNT 3

#endif

#ifdef SCAI_COMPLEX_SUPPORTED

typedef ComplexLongDouble ScalarRepType;

#else

typedef double ScalarRepType;

#endif

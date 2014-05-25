/**
 * @file LAMATypes.hpp
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
 * @brief LAMATypes.hpp
 * @author Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_LAMATYPES_HPP_
#define LAMA_LAMATYPES_HPP_

#include <cstring>
#include <limits>
#include <stdint.h>

/** LAMA uses for all its classes and routines an own namespace.
 *
 *  Applications using LAMA must either put \c lama:: before the used
 *  classes and routines or utilize the using directive of C++.
 *
 *  \code
 *                                         using namespace lama;
 *     lama::DenseVector<float> V;         DenseVector<float> V;
 *  \endcode
 */

namespace lama
{

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
const IndexType nIndex = std::numeric_limits<IndexType>::max();

/** Data type that is used for numbering of partitions.
 *
 */
typedef int PartitionId;

/** Definition for a constant value that indicates a non-available partition.
 */
const PartitionId nPartition = std::numeric_limits<PartitionId>::max();

/** Enumeration type for storage order of two-dimensional arrays, taken
 *  over from cblas.
 */

enum CBLAS_ORDER
{
    CblasRowMajor = 101, CblasColMajor = 102
};

/** Enumeration type for transpose use of two-dimensional arrays, taken
 *  over from cblas.
 */

enum CBLAS_TRANSPOSE
{
    CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113
};

/** Enumeration type for partial use of two-dimensional arrays, taken
 *  over from cblas.
 */
enum CBLAS_UPLO
{
    CblasUpper = 121, CblasLower = 122
};

enum CBLAS_DIAG
{
    CblasNonUnit = 131, CblasUnit = 132
};

enum CBLAS_SIDE
{
    CblasLeft = 141, CblasRight = 142
};

}

// Number of supported arithmetic types, maximal number is currently 4

#define ARITHMETIC_TYPE_CNT 3

// List here all arithmetic types for which matrices, storages might be created

#define ARITHMETIC_TYPE0        float
#define ARITHMETIC_TYPE1        double
#define ARITHMETIC_TYPE2        LongDouble

// Define for the arithmetic types the counterparts of enum Scalar::Tyep

#define SCALAR_ARITHMETIC_TYPE0 FLOAT
#define SCALAR_ARITHMETIC_TYPE1 DOUBLE
#define SCALAR_ARITHMETIC_TYPE2 LONG_DOUBLE

// For convenience we define ARRAY_TYPE

#define ARRAY_TYPE_CNT 4

#define ARRAY_TYPE0    int
#define ARRAY_TYPE1    float
#define ARRAY_TYPE2    double
#define ARRAY_TYPE3    LongDouble


#endif // LAMA_LAMATYPES_HPP_

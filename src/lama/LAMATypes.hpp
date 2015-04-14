/**
 * @file LAMATypes.hpp
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
 * @brief LAMATypes.hpp
 * @author Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_LAMATYPES_HPP_
#define LAMA_LAMATYPES_HPP_

// include LAMA
#include <lama/Complex.hpp>

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

/** Data type for complex numbers in single precision.
 *  LAMA uses its own data type instead of std::complex.
 */
typedef Complex<float> ComplexFloat;

/** Data type for complex numbers in double precision.
 *  LAMA uses its own data type instead of std::complex.
 */
typedef Complex<double> ComplexDouble;

/** Data type for complex numbers in long double precision.
 *  LAMA uses its own data type instead of std::complex.
 */
typedef Complex<long double> ComplexLongDouble;

/** Definition for a constant value that indicates a non-available index.
 */
static const IndexType nIndex = std::numeric_limits<IndexType>::max();

/** Data type that is used for numbering of partitions.
 *
 */
typedef int PartitionId;

/** Definition for a constant value that indicates a non-available partition.
 */
static const PartitionId nPartition = std::numeric_limits<PartitionId>::max();

} // namespace lama

// Number of supported arithmetic types, maximal number is currently 4

#define ARITHMETIC_TYPE_CNT 3

// List here all arithmetic types for which matrices, storages might be created

#define ARITHMETIC_TYPE0        float
#define ARITHMETIC_TYPE1        double
#define ARITHMETIC_TYPE2        ComplexFloat
#define ARITHMETIC_TYPE3        LongDouble

// Define for the arithmetic types the counterparts of enum Scalar::Tyep
// Sorry, we cannot use the routine getType<ARITHMETIC_TYPE##I> in case stmt

#define SCALAR_ARITHMETIC_TYPE0 FLOAT
#define SCALAR_ARITHMETIC_TYPE1 DOUBLE
#define SCALAR_ARITHMETIC_TYPE2 COMPLEX
#define SCALAR_ARITHMETIC_TYPE3 LONG_DOUBLE

// For convenience we define ARRAY_TYPE

#define ARRAY_TYPE_CNT 4

#define ARRAY_TYPE0    int
#define ARRAY_TYPE1    float
#define ARRAY_TYPE2    double
#define ARRAY_TYPE3    ComplexFloat
#define ARRAY_TYPE4    LongDouble

#endif // LAMA_LAMATYPES_HPP_

/**
 * @file SCAITypes.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Define all arithmetic types for which SCAI template classes will be instantiated
 * @author Jiri Kraus
 * @date 23.02.2011
 */
#pragma once

// local library

// no support of Complex if this file is not included

#ifdef SCAI_COMPLEX_SUPPORTED
    #include <scai/common/Complex.hpp>
#endif

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
/*
 * Lists of all used arithmetic types, two possibilities with/without complex
 */
#ifdef SCAI_COMPLEX_SUPPORTED
    /*
     * Kernel
     */
    #define SCAI_ARITHMETIC_HOST_CNT 6
    #define SCAI_ARITHMETIC_HOST     float, double, ComplexFloat, ComplexDouble, long double, ComplexLongDouble

    #define SCAI_ARITHMETIC_EXT_HOST_CNT 4
    #define SCAI_ARITHMETIC_EXT_HOST float, double, ComplexFloat, ComplexDouble

    #define SCAI_ARITHMETIC_CUDA_CNT 4
    #define SCAI_ARITHMETIC_CUDA     float, double, ComplexFloat, ComplexDouble

    #define SCAI_ARITHMETIC_MIC_CNT 4
    #define SCAI_ARITHMETIC_MIC      float, double, ComplexFloat, ComplexDouble

    /*
     * Array
     */

    #define SCAI_ARITHMETIC_ARRAY_HOST_CNT 7
    #define SCAI_ARITHMETIC_ARRAY_HOST     IndexType, SCAI_ARITHMETIC_HOST

    #define SCAI_ARITHMETIC_ARRAY_EXT_HOST_CNT 5
    #define SCAI_ARITHMETIC_ARRAY_EXT_HOST IndexType, SCAI_ARITHMETIC_EXT_HOST

    #define SCAI_ARITHMETIC_ARRAY_CUDA_CNT 5
    #define SCAI_ARITHMETIC_ARRAY_CUDA     IndexType, SCAI_ARITHMETIC_CUDA

    #define SCAI_ARITHMETIC_ARRAY_MIC_CNT 5
    #define SCAI_ARITHMETIC_ARRAY_MIC      IndexType, SCAI_ARITHMETIC_MIC
#else
    /*
     * Kernel
     */
    #define SCAI_ARITHMETIC_HOST_CNT 3
    #define SCAI_ARITHMETIC_HOST     float, double, long double

    #define SCAI_ARITHMETIC_EXT_HOST_CNT 2
    #define SCAI_ARITHMETIC_EXT_HOST float, double

    #define SCAI_ARITHMETIC_CUDA_CNT 2
    #define SCAI_ARITHMETIC_CUDA     float, double

    #define SCAI_ARITHMETIC_MIC_CNT 2
    #define SCAI_ARITHMETIC_MIC      float, double

    /*
     * Array
     */

    #define SCAI_ARITHMETIC_ARRAY_HOST_CNT 4
    #define SCAI_ARITHMETIC_ARRAY_HOST     IndexType, SCAI_ARITHMETIC_HOST

    #define SCAI_ARITHMETIC_ARRAY_EXT_HOST_CNT 3
    #define SCAI_ARITHMETIC_ARRAY_EXT_HOST IndexType, SCAI_ARITHMETIC_EXT_HOST

    #define SCAI_ARITHMETIC_ARRAY_CUDA_CNT 3
    #define SCAI_ARITHMETIC_ARRAY_CUDA     IndexType, SCAI_ARITHMETIC_CUDA

    #define SCAI_ARITHMETIC_ARRAY_MIC_CNT 3
    #define SCAI_ARITHMETIC_ARRAY_MIC      IndexType, SCAI_ARITHMETIC_MIC
#endif


/*
 * List creation
 */

#define SCAI_ARITHMETIC_HOST_LIST TYPELIST(     SCAI_ARITHMETIC_HOST_CNT,     SCAI_ARITHMETIC_HOST )
#define SCAI_ARITHMETIC_EXT_HOST_LIST TYPELIST( SCAI_ARITHMETIC_EXT_HOST_CNT, SCAI_ARITHMETIC_EXT_HOST )
#define SCAI_ARITHMETIC_CUDA_LIST TYPELIST(     SCAI_ARITHMETIC_CUDA_CNT,     SCAI_ARITHMETIC_CUDA )
#define SCAI_ARITHMETIC_MIC_LIST TYPELIST(      SCAI_ARITHMETIC_MIC_CNT,      SCAI_ARITHMETIC_MIC )

#define SCAI_ARITHMETIC_ARRAY_HOST_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_HOST_CNT, SCAI_ARITHMETIC_ARRAY_HOST )
#define SCAI_ARITHMETIC_ARRAY_EXT_HOST_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_EXT_HOST_CNT, SCAI_ARITHMETIC_ARRAY_EXT_HOST )
#define SCAI_ARITHMETIC_ARRAY_CUDA_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_CUDA_CNT, SCAI_ARITHMETIC_ARRAY_CUDA )
#define SCAI_ARITHMETIC_ARRAY_MIC_LIST TYPELIST( SCAI_ARITHMETIC_ARRAY_MIC_CNT, SCAI_ARITHMETIC_ARRAY_MIC )

/** Number of supported types used in REPEAT macros */

#ifdef SCAI_COMPLEX_SUPPORTED

typedef ComplexLongDouble ScalarRepType;

#else

typedef double ScalarRepType;

#endif

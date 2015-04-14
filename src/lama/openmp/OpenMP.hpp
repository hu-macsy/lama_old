/**
 * @file OpenMP.hpp
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
 * @brief Common defintions for optional use of OpenMP
 * @author Thomas Brandes
 * @date 11.06.2013
 * @since 1.0.1
 */
#ifndef LAMA_OPENMP_HPP_
#define LAMA_OPENMP_HPP_

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
    #define omp_get_max_threads() 1
#endif

/** atomicAdd used for reductions as reduction directive is unsupported for complex numbers.
 *
 *  Note: template specialization used for float and double
 */

template<typename ValueType>
inline void atomicAdd( ValueType& sharedResult, const ValueType& threadResult )
{
#pragma omp critical
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( float& sharedResult, const float& threadResult )
{
#pragma omp atomic
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( double& sharedResult, const double& threadResult )
{
#pragma omp atomic
    sharedResult += threadResult;
}

#endif //  LAMA_OPENMP_HPP_

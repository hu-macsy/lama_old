/**
 * @file PMVBenchmark.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief PMVBenchmark.cpp
 * @author brandes
 * @date 11.05.2011
 * $Id$
 */

#include <bench/PMVBenchmark.hpp>

#include <framework/src/benchmark_framework.hpp>

// Matrix storage types needed for instantiation

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#define LAMAMVBENCHMARKREGISTRATION( MATRIX, GROUP )                            \
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    PMVBenchmark<MATRIX<float> >::group()                                           \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                 \
        return group;                                                               \
    }                                                                               \
    \
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    PMVBenchmark<MATRIX<double> >::group()                                          \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                 \
        return group;                                                               \
    }                                                                               \
    \
    template<>                                                                      \
    const std::string& PMVBenchmark<MATRIX<float> >::sid()                          \
    {                                                                               \
        static const std::string sid = "LAMA<float>";                               \
        return sid;                                                                 \
    }                                                                               \
    \
    template<>                                                                      \
    const std::string& PMVBenchmark<MATRIX<double> >::sid()                         \
    {                                                                               \
        static const std::string sid = "LAMA<double>";                              \
        return sid;                                                                 \
    }                                                                               \
    \
    \
    LAMA_BENCHMARK_REGISTRATION2(PMVBenchmark<MATRIX<float> >,GROUP##float);             \
    LAMA_BENCHMARK_REGISTRATION2(PMVBenchmark<MATRIX<double> >,GROUP##double);

LAMAMVBENCHMARKREGISTRATION( lama::CSRSparseMatrix, CSRSAMGSpMV )
LAMAMVBENCHMARKREGISTRATION( lama::ELLSparseMatrix, ELLSpMV )
LAMAMVBENCHMARKREGISTRATION( lama::COOSparseMatrix, COOSpMV )
LAMAMVBENCHMARKREGISTRATION( lama::DIASparseMatrix, DIASpMV )
LAMAMVBENCHMARKREGISTRATION( lama::JDSSparseMatrix, JDSSpMV )
LAMAMVBENCHMARKREGISTRATION( lama::DenseMatrix, DenseMV )

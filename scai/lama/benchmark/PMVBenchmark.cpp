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

#include <scai/lama/benchmark/PMVBenchmark.hpp>

// Matrix storage types needed for instantiation

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#define LAMAMVBENCHMARKREGISTRATION( MATRIX, GROUP )                                \
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    PMVBenchmark<MATRIX<float> >::group()                                           \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                     \
        return group;                                                               \
    }                                                                               \
                                                                                    \
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    PMVBenchmark<MATRIX<double> >::group()                                          \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                     \
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
    LAMA_BENCHMARK_REGISTRATION2(PMVBenchmark<MATRIX<float> >,GROUP##float);        \
    LAMA_BENCHMARK_REGISTRATION2(PMVBenchmark<MATRIX<double> >,GROUP##double);

LAMAMVBENCHMARKREGISTRATION( scai::lama::CSRSparseMatrix, CSRSAMGSpMV )
LAMAMVBENCHMARKREGISTRATION( scai::lama::ELLSparseMatrix, ELLSpMV )
LAMAMVBENCHMARKREGISTRATION( scai::lama::COOSparseMatrix, COOSpMV )
LAMAMVBENCHMARKREGISTRATION( scai::lama::DIASparseMatrix, DIASpMV )
LAMAMVBENCHMARKREGISTRATION( scai::lama::JDSSparseMatrix, JDSSpMV )
LAMAMVBENCHMARKREGISTRATION( scai::lama::DenseMatrix, DenseMV )

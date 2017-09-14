/**
 * @file PMVBenchmark.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief PMVBenchmark.cpp
 * @author brandes
 * @date 11.05.2011
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

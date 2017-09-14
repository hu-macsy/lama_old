/**
 * @file ConvertMatrixStorageBenchmark.cpp
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
 * @brief ConvertMatrixStorageBenchmark.cpp
 * @author Jiri Kraus and Bea Hornef
 * @date 02.12.2011
 */

#include <scai/lama/benchmark/ConvertMatrixStorageBenchmark.hpp>

#include <scai/lama/benchmark/LAMAInputSetComplexityVisitor.hpp>

// Matrix storage types needed for instantiation

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>

#define LAMACONVERTMATRIXSTORAGEBENCHMARKREGISTRATION( STORAGE1, STORAGE2, GROUP )               \
    template<>                                                                                   \
    const LAMAInputSetComplexityVisitor::Group&                                                  \
    ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> >::group()                    \
    {                                                                                            \
        static const LAMAInputSetComplexityVisitor::Group group                                  \
        = LAMAInputSetComplexityVisitor::GROUP;                                                  \
        return group;                                                                            \
    }                                                                                            \
    \
    template<>                                                                                   \
    const LAMAInputSetComplexityVisitor::Group&                                                  \
    ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> >::group()                  \
    {                                                                                            \
        static const LAMAInputSetComplexityVisitor::Group group                                  \
        = LAMAInputSetComplexityVisitor::GROUP;                                                  \
        return group;                                                                            \
    }                                                                                            \
    \
    template<>                                                                                   \
    const std::string& ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> >::sid()   \
    {                                                                                            \
        static const std::string sid = "LAMA<float>";                                            \
        return sid;                                                                              \
    }                                                                                            \
    \
    template<>                                                                                   \
    const std::string& ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> >::sid() \
    {                                                                                            \
        static const std::string sid = "LAMA<double>";                                           \
        return sid;                                                                              \
    }                                                                                            \
    \
    typedef ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> > MYTYPE##float;      \
    typedef ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> > MYTYPE##double;   \
    \
    LAMA_BENCHMARK_REGISTRATION2(MYTYPE##float,GROUP##float);                                    \
    LAMA_BENCHMARK_REGISTRATION2(MYTYPE##double,GROUP##double);

LAMACONVERTMATRIXSTORAGEBENCHMARKREGISTRATION( lama::CSRStorage, lama::JDSStorage, ConvertCSR2JDS )


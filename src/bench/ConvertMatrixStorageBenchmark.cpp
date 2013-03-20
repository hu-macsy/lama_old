/**
 * @file ConvertMatrixStorageBenchmark.cpp
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
 * @brief ConvertMatrixStorageBenchmark.cpp
 * @author Jiri Kraus and Bea Hornef
 * @date 02.12.2011
 * $Id$
 */

#include <bench/ConvertMatrixStorageBenchmark.hpp>

#include <framework/src/benchmark_framework.hpp>

#include <bench/LAMAInputSetComplexityVisitor.hpp>
// Matrix storage types needed for instantiation

#include <lama/storage/CSRStorage.hpp>
#include <lama/storage/JDSStorage.hpp>

#define LAMACONVERTMATRIXSTORAGEBENCHMARKREGISTRATION( STORAGE1, STORAGE2, GROUP )\
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> >::group()       \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                 \
        return group;                                                               \
    }                                                                               \
    \
    template<>                                                                      \
    const LAMAInputSetComplexityVisitor::Group&                                     \
    ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> >::group()     \
    {                                                                               \
        static const LAMAInputSetComplexityVisitor::Group group                     \
        = LAMAInputSetComplexityVisitor::GROUP;                                 \
        return group;                                                               \
    }                                                                               \
    \
    template<>                                                                      \
    const std::string& ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> >::sid()\
    {                                                                               \
        static const std::string sid = "LAMA<float>";                               \
        return sid;                                                                 \
    }                                                                               \
    \
    template<>                                                                      \
    const std::string& ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> >::sid()\
    {                                                                               \
        static const std::string sid = "LAMA<double>";                              \
        return sid;                                                                 \
    }                                                                               \
    \
    typedef ConvertMatrixStorageBenchmark<STORAGE1<float>, STORAGE2<float> > MYTYPE##float;\
    typedef ConvertMatrixStorageBenchmark<STORAGE1<double>, STORAGE2<double> > MYTYPE##double;\
    \
    LAMA_BENCHMARK_REGISTRATION2(MYTYPE##float,GROUP##float);            \
    LAMA_BENCHMARK_REGISTRATION2(MYTYPE##double,GROUP##double);

LAMACONVERTMATRIXSTORAGEBENCHMARKREGISTRATION( lama::CSRStorage, lama::JDSStorage, ConvertCSR2JDS )


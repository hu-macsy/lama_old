/**
 * @file kregistry/examples/BenchInterface.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Benchmark that measures the advantage of static variables for kernel functions.
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>

using namespace scai;
using namespace scai::kregistry;

using scai::common::ContextType;

template<typename ValueType>
static ValueType add( ValueType x )
{
    return x + 1;
}

template<typename ValueType>
static ValueType sub( ValueType x )
{
    return x - 1;
}

const char* add_names[] = { "A+", "B+", "C+", "D+", "E+", "F+", "G+", "H+", "I+", "J+",
                            "K+", "L+", "M+", "N+", "O+", "P+", "Q+", "R+", "S+", "T+"
                          };

const char* sub_names[] = { "A-", "B-", "C-", "D-", "E-", "F-", "G-", "H-", "I-", "J-",
                            "K-", "L-", "M-", "N-", "O-", "P-", "Q-", "R-", "S-", "T-"
                          };

static void setInterface()
{
    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;

    // register 20  x 4 routines at the kernel registry

    for ( int i = 0; i < 20; ++i )
    {
        KernelRegistry::set( add<float>, add_names[i], ContextType::Host, flag );
        KernelRegistry::set( add<double>, add_names[i], ContextType::Host, flag );
        KernelRegistry::set( sub<float>, sub_names[i], ContextType::Host, flag );
        KernelRegistry::set( sub<double>, sub_names[i], ContextType::Host, flag );
    }
}

static void doIt1 ( double x )
{
    // Usual declaration, the functions are searched with each call
    KernelContextFunction< double (* ) ( double ) > add( "E+" );
    KernelContextFunction< double (* ) ( double ) > sub( "S-" );
    x = add[ContextType::Host]( sub[ContextType::Host]( x ) );
}

static void doIt2 ( double x )
{
    // static declaration, the functions are searched only in first call
    static KernelContextFunction< double (* ) ( double ) > add( "E+" );
    static KernelContextFunction< double (* ) ( double ) > sub( "S-" );
    x = add[ContextType::Host]( sub[ContextType::Host]( x ) );
}

int main()
{
    using scai::common::Walltime;
    setInterface();
    KernelRegistry::printAll();
    double x = 0.0;
    const int N = 100 * 1000;
    // measure for routine where kernel functions are search each time
    double time1 = Walltime::get();

    for ( int i = 0; i < N; ++ i )
    {
        doIt1( x );
    }

    time1 = Walltime::get() - time1;
    std::cout << "time1 = " << time1 * 1000.0 << " ms " << std::endl;
    SCAI_ASSERT_EQUAL( 0, x, "Wrong result" )
    double time2 = Walltime::get();

    for ( int i = 0; i < N; ++ i )
    {
        doIt2( x );
    }

    // measure for routine where kernel functions are looked up only once
    time2 = Walltime::get() - time2;
    std::cout << "time2 = " << time2 * 1000.0 << " ms " << std::endl;
    SCAI_ASSERT_EQUAL( 0, x, "Wrong result" )
    std::cout << "final x = " << x << ", should be 0.0" << std::endl;
    double c1_us = time1 * 1000.0 * 1000.0 / N;
    double c2_us = time2 * 1000.0 * 1000.0 / N;
    std::cout << "Summary ( N = " << N << " ) : dyn : " << c1_us << " us, stat: " << c2_us
              << ", ratio = " << ( c1_us / c2_us ) << std::endl;
}

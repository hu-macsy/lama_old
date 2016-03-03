/**
 * @file kregistry/examples/BenchInterface.cpp
 *
 * @brief Benchmark that measures the advantage of static variables for kernel functions.
 *
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>

using namespace scai;
using namespace scai::kregistry;

using scai::common::context;

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
                            "K+", "L+", "M+", "N+", "O+", "P+", "Q+", "R+", "S+", "T+" };

const char* sub_names[] = { "A-", "B-", "C-", "D-", "E-", "F-", "G-", "H-", "I-", "J-", 
                            "K-", "L-", "M-", "N-", "O-", "P-", "Q-", "R-", "S-", "T-" };

static void setInterface()
{
    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;

    // register 20  x 4 routines at the kernel registry

    for ( int i = 0; i < 20; ++i )
    {
        KernelRegistry::set( add<float>, add_names[i], context::Host, flag );
        KernelRegistry::set( add<double>, add_names[i], context::Host, flag );
        KernelRegistry::set( sub<float>, sub_names[i], context::Host, flag );
        KernelRegistry::set( sub<double>, sub_names[i], context::Host, flag );
    }
}

static void doIt1 ( double x ) 
{
    // Usual declaration, the functions are searched with each call

    KernelContextFunction< double (*) ( double ) > add( "E+" );
    KernelContextFunction< double (*) ( double ) > sub( "S-" );

    x = add[context::Host]( sub[context::Host]( x ) );
}

static void doIt2 ( double x ) 
{
    // static declaration, the functions are searched only in first call

    static KernelContextFunction< double (*) ( double ) > add( "E+" );
    static KernelContextFunction< double (*) ( double ) > sub( "S-" );

    x = add[context::Host]( sub[context::Host]( x ) );
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

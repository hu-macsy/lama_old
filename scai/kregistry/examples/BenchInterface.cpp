
#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

using namespace scai;
using namespace scai::kregistry;

using namespace scai::common::context;

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
 
const char* names[] = { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
                        "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T" };

static void setInterface()
{
    for ( int i = 0; i < 20; ++i )
    {
        KernelRegistry::set(add<float>, names[i], Host );
        KernelRegistry::set(add<double>, names[i], Host );
        KernelRegistry::set(sub<float>, names[i], Host );
        KernelRegistry::set(sub<double>, names[i], Host );
    }
}

static void doIt1 ( double x ) 
{
    KernelContextFunction< double (*) ( double ) > add( "E" );
    KernelContextFunction< double (*) ( double ) > sub( "S" );

    x = add[Host]( sub[Host]( x ) );
}

static void doIt2 ( double x ) 
{
    static KernelContextFunction< double (*) ( double ) > add( "E" );
    static KernelContextFunction< double (*) ( double ) > sub( "S" );

    x = add[Host]( sub[Host]( x ) );
}

int main()
{
    using scai::common::Walltime;

    setInterface();

    KernelRegistry::printAll();
 
    double x = 0.0;

    const int N = 100 * 1000;

    double time1 = Walltime::get();

    for ( int i = 0; i < N; ++ i )
    {
        doIt1( x );
    }

    time1 = Walltime::get() - time1;
  
    std::cout << "time1 = " << time1 * 1000.0 << " ms " << std::endl;

    double time2 = Walltime::get();

    for ( int i = 0; i < N; ++ i )
    {
        doIt2( x );
    }

    time2 = Walltime::get() - time2;
  
    std::cout << "time2 = " << time2 * 1000.0 << " ms " << std::endl;
    
    std::cout << "final x = " << x << ", should be 0.0" << std::endl;

    double c1_us = time1 * 1000.0 * 1000.0 / N;
    double c2_us = time2 * 1000.0 * 1000.0 / N;

    std::cout << "Summary ( N = " << N << " ) : dyn : " << c1_us << " us, stat: " << c2_us 
              << ", ratio = " << ( c1_us / c2_us ) << std::endl;
}

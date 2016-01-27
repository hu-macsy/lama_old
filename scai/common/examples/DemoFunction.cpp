/**
 * @file common/examples/DemoFunction.cpp
 *
 * @brief Example for using function / bind
 *
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/function.hpp>
#include <scai/common/bind.hpp>

int f ( double x, int p )
{
   double f = 1.0;

   for ( int i = 0; i < p; ++i )
   {
      f *= 10.0;
   }
   return static_cast<int>( x * f );
}

using namespace scai::common;

int main()
{
    int ( *foo1 ) ( double, int ) = &f;       // traditional function pointer
    function<int( double, int )> foo2 = &f;    // function wrapper class

    std::cout << "Call foo1( 3.15143, 2 ) = " << foo1( 3.15143, 2 ) << std::endl;
    std::cout << "Call foo2( 3.15143, 3 ) = " << foo2( 3.15143, 3 ) << std::endl;

    function<int( double )> foo2a = bind( f, _1, 5 );

    std::cout << "Call foo2a( 3.15143 ) = " << foo2a( 3.15143 )  << std::endl;

    function<int( int )> foo2b = bind( f, 3.15143, _1 );

    std::cout << "Call foo2b( 1 ) = " << foo2b( 1 )  << std::endl;

    function<int( int, double )> foo3 = bind( f, _2, _1 );
    std::cout << "Call foo3( 2, 3.15143 ) = " << foo3( 2, 3.15143 ) << std::endl;
}

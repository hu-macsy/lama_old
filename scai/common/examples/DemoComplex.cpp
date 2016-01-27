/**
 * @file common/examples/DemoComplex.cpp
 *
 * @brief Using complex numbers of scai/common
 *
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/Complex.hpp>
#include <complex>

#include <iostream>

using scai::common::Complex;
using scai::common::conj;

int main()
{
    Complex<float> a( 1.0f, 1.0f );
    Complex<double> b( 2.0, -1.0 );
    std::cout << "a = " << a << ", b = " << b << std::endl;
 
    double b_re = b.real();
    double b_im = b.imag();

    std::cout << "b.real() = " << b_re << ", b.imag() = " << b_im << std::endl;
    std::cout << "conj( b ) = " << conj( b ) << std::endl;
    // not allowed for std::complex
    std::cout << "static_cast<double>( b ) = " << static_cast<double>( b ) << std::endl;
    std::cout << "abs( b ) = " << abs( b ) << std::endl;

    if ( a < b )
    {
        std::cout << a << " < " << b << std::endl;
    }
    else if ( a == b )
    {
        std::cout << a << " == " << b << std::endl;
    }
    else if ( a >= b )
    {
        std::cout << a << " >= " << b << std::endl;
    }
}

/**
 * @file common/examples/DemoTypeTrait.cpp
 *
 * @brief Examples for using TypeTraits
 *
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/Complex.hpp>
#include <scai/common/Math.hpp>

#include <iostream>

using scai::common::Math;

template<typename ValueType>
void testRoutine( ValueType x )
{
    ValueType absX = Math::abs( x );

    std::cout << "abs ( " << x << " ) = " << absX << std::endl;

    ValueType sqrtX = Math::sqrt( x );

    std::cout << "sqrt ( " << x << " ) = " << sqrtX << std::endl;

    ValueType conjX = Math::conj( x );

    std::cout << "conj ( " << x << " ) = " << conjX << std::endl;
}

int main()
{
    testRoutine<float>( -1 );
    testRoutine<double>( -1 );
    testRoutine<ComplexFloat>( ComplexFloat( 2, -1 ) );
}

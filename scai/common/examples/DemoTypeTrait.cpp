/**
 * @file common/examples/DemoTypeTrait.cpp
 *
 * @brief Examples for using TypeTraits
 *
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/TypeTraits.hpp>

using scai::common::TypeTraits;

template<typename ValueType>
void testRoutine()
{
    std::cout << "TypeTraits<...>::id() = " << TypeTraits<ValueType>::id() << std::endl;

    //  problem here
    //  std::cout << "TypeTraits<...>::stype = " << TypeTraits<ValueType>::stype << std::endl;

    std::cout << "TypeTraits<...>::stype = " << scai::common::getScalarType<ValueType>() << std::endl;

    ValueType x( -1 );

    ValueType absX = TypeTraits<ValueType>::abs( x );

    std::cout << "abs ( " << x << " ) = " << absX << std::endl;

    ValueType sqrtX = TypeTraits<ValueType>::sqrt( x );

    std::cout << "sqrt ( " << x << " ) = " << sqrtX << std::endl;
}

int main()
{
    testRoutine<float>();
    testRoutine<double>();
    testRoutine<ComplexFloat>();
}

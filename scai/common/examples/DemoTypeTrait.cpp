/**
 * @file common/examples/DemoTypeTrait.cpp
 *
 * @brief Examples for using TypeTraits
 *
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>

using scai::common::TypeTraits;

template<typename ValueType>
void testRoutine()
{
    std::cout << "TypeTraits<...>::id() = " << TypeTraits<ValueType>::id() << std::endl;

    scai::common::scalar::ScalarType x = TypeTraits<ValueType>::stype;
    std::cout << "TypeTraits<...>::stype = " << x << std::endl;
}

int main()
{
    testRoutine<float>();
    testRoutine<double>();
    testRoutine<ComplexFloat>();
}

/**
 * @file common/examples/DemoTypeTrait.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Examples for using TypeTraits
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>

#include<iostream>
#include<iomanip>

using scai::common::TypeTraits;

template<typename ValueType>
void testRoutine()
{
    std::cout << "TypeTraits<...>::id() = " << TypeTraits<ValueType>::id() << std::endl;
    scai::common::scalar::ScalarType stype = TypeTraits<ValueType>::stype;
    std::cout << "TypeTraits<...>::stype = " << stype << std::endl;
    std::cout << "isComplex = " << isComplex( stype ) << std::endl;

    ValueType alpha = ValueType( 1 ) / ValueType( 3 );
    int precision = TypeTraits<ValueType>::precision();
    std::cout << "Output alpha: precision = " << precision;
    std::cout << std::setprecision( precision ) << " alpha = " << alpha << std::endl;
  
}

int main()
{
    testRoutine<float>();
    testRoutine<double>();
    testRoutine<long double>();
#ifdef SCAI_COMPLEX_SUPPORTED
    testRoutine<ComplexFloat>();
#endif
}

/**
 * @file common/examples/DemoTypeTrait.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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

using namespace scai;
using common::TypeTraits;

using std::cout;
using std::endl;
using std::string;
using std::setprecision;

template<typename ValueType>
void testRoutine()
{
    // common string used for output of name of ValueType

    string typeTemplate = TypeTraits<ValueType>::id();
    typeTemplate = "<" + typeTemplate + ">";

    cout << "TypeTraits" << typeTemplate << " = " << TypeTraits<ValueType>::id() << endl;
    common::ScalarType stype = TypeTraits<ValueType>::stype;
    cout << "TypeTraits" << typeTemplate << "::stype = " << stype << endl;
    cout << "isComplex = " << isComplex( stype ) << endl;

    ValueType alpha = ValueType( 1 ) / ValueType( 3 );
    ValueType beta  = ValueType( 2 ) / ValueType( 3 );
    int precision = TypeTraits<ValueType>::precision();
    cout << "precision" << typeTemplate << " = " << precision << endl;
    cout << setprecision( precision ) << "  1/3 = " << alpha << endl;
    cout << setprecision( precision ) << "  2/3 = " << beta  << endl;

    cout << "eps0" << typeTemplate << " = " << TypeTraits<ValueType>::eps0() << endl;
    cout << "eps1" << typeTemplate << " = " << TypeTraits<ValueType>::eps1() << endl;
    cout << "min" << typeTemplate << " = " << TypeTraits<ValueType>::getMin() << endl;
    cout << "max" << typeTemplate << " = " << TypeTraits<ValueType>::getMax() << endl;
    cout << "small" << typeTemplate << " = " << TypeTraits<ValueType>::small() << endl;
}

int main()
{
    testRoutine<IndexType>();
    testRoutine<float>();
    testRoutine<double>();
    testRoutine<long double>();

#ifdef SCAI_COMPLEX_SUPPORTED
    testRoutine<ComplexFloat>();
    testRoutine<ComplexDouble>();
    testRoutine<ComplexLongDouble>();
#endif

}

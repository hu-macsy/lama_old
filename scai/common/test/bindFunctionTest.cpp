/**
 * @file bindFunctionTest.cpp
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
 * @brief Test routines for bind and function wrapper
 * @author Lauretta Schubert
 * @date 30.03.2016
 */

// follows the c++ example (http://www.cplusplus.com/reference/functional/bind/)

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/function.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace common;

/* -------------------------------------------------------------------------------- */

template<typename T>
T my_divide (T x, T y)
{
    return x / y;
}

template<typename T>
struct MyPair
{
    MyPair( T _a, T _b ): a(_a), b(_b) {}
    T a,b;
    T multiply()
    {
        return a * b;
    }
};

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( bindFunctionTest, ValueType, scai_arithmetic_test_types )
{
    function<ValueType()> fn_five = bind ( my_divide<ValueType>, static_cast<ValueType>(10),static_cast<ValueType>(2) );
    ValueType res1 = fn_five();

    BOOST_CHECK_EQUAL( res1, 5.0 );

    MyPair<ValueType> ten_two( 10, 2 );
    function<ValueType( MyPair<ValueType> )> bound_member_fn = bind ( &MyPair<ValueType>::multiply, _1 );
    ValueType res2 = bound_member_fn( ten_two );

    BOOST_CHECK_EQUAL( res2, 20.0 );
}

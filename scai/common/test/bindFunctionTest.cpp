/**
 * @file bindFunctionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Test routines for bind and function wrapper
 *
 * @author Lauretta Schubert
 * @date 30.03.2016
 */

// follows the c++ example (http://www.cplusplus.com/reference/functional/bind/)

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/function.hpp>

using namespace scai;
using namespace common;

/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<ARITHMETIC_HOST> SCAI_ARITHMETIC_TYPES;

/* -------------------------------------------------------------------------------- */

template<typename T>
T my_divide (T x, T y) { return x / y; }

template<typename T>
struct MyPair {
  MyPair( T _a, T _b ): a(_a), b(_b) {}
  T a,b;
  T multiply() { return a * b; }
};

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( bindFunctionTest, ValueType, SCAI_ARITHMETIC_TYPES )
{
	function<ValueType()> fn_five = bind ( my_divide<ValueType>, static_cast<ValueType>(10),static_cast<ValueType>(2) );
	ValueType res1 = fn_five();

	BOOST_CHECK_EQUAL( res1, 5.0 );

	MyPair<ValueType> ten_two( 10, 2 );
	function<ValueType( MyPair<ValueType> )> bound_member_fn = bind ( &MyPair<ValueType>::multiply, _1 );
	ValueType res2 = bound_member_fn( ten_two );

	BOOST_CHECK_EQUAL( res2, 20.0 );
}

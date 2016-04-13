/**
 * @file shared_ptrTest.cpp
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
 * @brief Test routines for smart pointer wrapper
 *
 * @author Lauretta Schubert
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/weak_ptr.hpp>

using namespace scai;
using namespace common;

/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_HOST> ValueTypes;

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( shared_ptrTest, ValueType, ValueTypes )
{
	shared_ptr<ValueType> emptyPointer;

	BOOST_CHECK_EQUAL ( emptyPointer.use_count(), 0 );

  	shared_ptr<ValueType> pointer ( new ValueType(10) );

    BOOST_CHECK_EQUAL( pointer.use_count(), 1 );

    shared_ptr<ValueType> secPointer ( pointer );

    BOOST_CHECK_EQUAL( secPointer.use_count(), 2 );
}

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unique_ptrTest, ValueType, ValueTypes )
{
	unique_ptr<ValueType> emptyPointer;

	bool test = ( emptyPointer.get() == NULL );

	BOOST_CHECK ( test );

	unique_ptr<ValueType> pointer ( new ValueType(10) );

	test = ( pointer.get() != NULL );

	BOOST_CHECK ( test );
}

/* -------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( weak_ptrTest, ValueType, ValueTypes )
{
	weak_ptr<ValueType> emptyPointer;

	BOOST_CHECK_EQUAL ( emptyPointer.use_count(), 0 );

	weak_ptr<ValueType> sec_emptyPointer( emptyPointer );

	BOOST_CHECK_EQUAL ( sec_emptyPointer.use_count(), 0 );	

  	shared_ptr<ValueType> shared_pointer ( new ValueType(10) );
    weak_ptr<ValueType> pointer ( shared_pointer );

    BOOST_CHECK_EQUAL( pointer.use_count(), 1 );

    shared_pointer.reset();

    BOOST_CHECK_EQUAL( pointer.use_count(), 0 );
}

/* -------------------------------------------------------------------------------- */

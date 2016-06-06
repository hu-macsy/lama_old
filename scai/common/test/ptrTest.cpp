/**
 * @file test/ptrTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Test routines for smart pointer wrapper
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

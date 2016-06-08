/**
 * @file TypeTraitsTest.cpp
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
 * @brief Test routines for class TypeTraits
 * @author Thomas Brandes
 * @date 05.02.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace common;

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( TypeTraitsTest, ValueType, scai_arithmetic_test_types )
{
    scalar::ScalarType stype = TypeTraits<ValueType>::stype;
    std::ostringstream out1;
    std::ostringstream out2;
    out1 << TypeTraits<ValueType>::id();
    out2 << stype;
    // std::cout << "out1 = " << out1.str() << ", out2 = " << out2.str() << std::endl;
    BOOST_CHECK_EQUAL( out1.str(), out2.str() );
}

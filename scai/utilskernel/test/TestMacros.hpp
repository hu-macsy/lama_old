/**
 * @file utilskernel/test/TestMacros.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Own SCAI test macros used for utilskernel
 * @author Thomas Brandes
 * @date 12.07.2018
 */

#include <scai/common/test/TestMacros.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

/** This macros checks that all elements of two arrays are the same  */

#define SCAI_CHECK_EQUAL_ARRAY( array1, array2 )                                                \
{                                                                                               \
    BOOST_TEST( scai::hmemo::hostReadAccess( array1 ) == scai::hmemo::hostReadAccess( array2 ), \
                boost::test_tools::per_element() );                                             \
}  

#define SCAI_CHECK_SMALL_ARRAY_DIFF( array1, array2, eps )                                          \
{                                                                                                   \
    BOOST_REQUIRE_EQUAL( array1.size(), array2.size() );                                            \
                                                                                                    \
    auto diff = scai::utilskernel::HArrayUtils::maxDiffNorm( array1, array2 );                      \
                                                                                                    \
    if ( diff > eps )                                                                               \
    {                                                                                               \
        BOOST_TEST( scai::hmemo::hostReadAccess( array1 ) == scai::hmemo::hostReadAccess( array2 ), \
                    boost::test_tools::per_element() );                                             \
    }                                                                                               \
}  

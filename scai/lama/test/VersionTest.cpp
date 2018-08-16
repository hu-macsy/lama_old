/**
 * @file VersionTest.cpp
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
 * @brief Contains test for Version check.
 * @author Thomas Brandes
 * @date 03.06.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama.hpp>
#include <string>

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VersionTest )
{
    std::string version = lama_get_version();
    // version is something like x.y.z
    int len = version.length();
    BOOST_CHECK( len >= 5 );
}


/**
 * @file lama_cuda_test.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class lama_cuda_test.
 * @author Alexander Büchel
 * @date 31.01.2012
 */
#ifndef WIN32
#define BOOST_ALL_DYN_LINK
#endif //WIN32
#define BOOST_TEST_MODULE lama_cuda_test

#include <boost/test/unit_test.hpp>

struct GlobalFixture
{
    GlobalFixture()
    {
    }
    ~GlobalFixture()
    {
    }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );


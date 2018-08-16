/**
 * @file scai/hmemo/test/ContextFix.hpp
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
 * @brief Test fixture for context
 * @author Thomas Brandes
 * @date 15.03.2016
 */

#pragma once

#include <scai/hmemo/Context.hpp>

#include <boost/test/unit_test.hpp>

/* --------------------------------------------------------------------- */

/** Fixture to be used for BOOST_GLOBAL_FIXTURE
 *
 *  provides access to testContext used as context at which tests should run
 *
 *  Purpose: use global Fixture avoids init/free of context device for each single test
 *
 *  Note: static variable ContextFix::testContext must be defined in cpp file.
 */
struct ContextFix
{
    ContextFix()
    {
        testContext = scai::hmemo::Context::getContextPtr();
        // BOOST_TEST_MESSAGE( "Setup ContextFix: test context = " << *testContext );
    }

    ~ContextFix()
    {
        // BOOST_TEST_MESSAGE( "Teardown ContextFix" );
        testContext.reset();
    }

    static scai::hmemo::ContextPtr testContext;
};


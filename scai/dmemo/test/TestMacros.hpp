/**
 * @file scai/dmemo/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <scai/common/test/TestMacros.hpp>

#define CHECK_COMMUNICATION_PLANS_EQUAL(plan1, plan2)                       \
    BOOST_TEST_CONTEXT(" CommunicationPlan instances do not match ")        \
    {                                                                       \
        BOOST_TEST(plan1.totalQuantity() == plan2.totalQuantity());         \
        BOOST_TEST(plan1.size() == plan2.size());                           \
                                                                            \
        if (plan1.size() == plan2.size())                                   \
        {                                                                   \
            for (PartitionId i = 0; i < plan1.size(); ++i)                  \
            {                                                               \
                BOOST_TEST_CONTEXT("mismatch at entry " << i)               \
                {                                                           \
                    const auto entry1 = plan1[i];                           \
                    const auto entry2 = plan2[i];                           \
                                                                            \
                    BOOST_TEST(entry1.partitionId == entry2.partitionId);   \
                    BOOST_TEST(entry1.quantity == entry2.quantity);         \
                    BOOST_TEST(entry1.offset == entry2.offset);             \
                }                                                           \
            }                                                               \
        }                                                                   \
        else                                                                \
        {                                                                   \
            BOOST_TEST(plan1.size() == plan2.size());                       \
        }                                                                   \
    }


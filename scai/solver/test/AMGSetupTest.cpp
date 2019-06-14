/**
 * @file AMGSetupTest.cpp
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
 * @brief Test of all available AMG setups from factory.
 * @author Thomas Brandes
 * @date 16.05.2019
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/AMGSetup.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <scai/hmemo/Context.hpp>

#include <scai/solver/test/TestMacros.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/LibModule.hpp>

using namespace scai;

using namespace lama;
using namespace solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( AMGSetupTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.AMGSetupTest" )

template<typename ValueType>
std::vector<std::shared_ptr<AMGSetup<ValueType>>> getAllAMGSetups()
{
    static bool initialized = false;  // make sure that dynamic libs will be loaded only once

    if ( !initialized )
    {
        std::string loadPath;

        if ( scai::common::Settings::getEnvironment( loadPath, "SCAI_LIBRARY_PATH" ) )
        {
            SCAI_LOG_INFO( logger, "Load all module libraries in loadPath " << loadPath )
            scai::common::LibModule::loadLibsByPath( loadPath.c_str() );
        }
        else
        {
            SCAI_LOG_WARN( logger, "No additional setup modules loaded, SCAI_LIBRARY_PATH not set" )
        }

        initialized = true;
    }

    std::vector<std::shared_ptr<AMGSetup<ValueType>>> amgSetups;

    std::vector<std::string> values;  // string is create type for the factory

    AMGSetup<ValueType>::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::shared_ptr<AMGSetup<ValueType> > setup( AMGSetup<ValueType>::getAMGSetup( values[i] ) );
        amgSetups.push_back( setup );
    }

    return amgSetups;
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    typedef DefaultReal ValueType;

    auto allSetups = getAllAMGSetups<ValueType>();

    for ( size_t i = 0; i < allSetups.size(); ++i )
    {
        auto& setup = *allSetups[i];

        SCAI_LOG_INFO( logger, "Test for AMG setup (uninitialized): " << setup )

        BOOST_CHECK_EQUAL( setup.getNumLevels(), 0 );

        CSRSparseMatrix<ValueType> inputA;

        const IndexType N = 25;

        MatrixCreator::buildPoisson2D( inputA, 9, N, N );

        setup.initialize( inputA );

        SCAI_LOG_INFO( logger, "Test for AMG setup (initialized): " << setup )

        IndexType numLevels = setup.getNumLevels();

        BOOST_CHECK( numLevels > 0 );

        for ( IndexType i = 0; i < numLevels; ++i )
        {
            const Matrix<ValueType>& galerkin = setup.getGalerkin( i );
            BOOST_CHECK_EQUAL( galerkin.getRowDistribution(), galerkin.getColDistribution() );
            BOOST_CHECK( galerkin.getNumRows() <= inputA.getNumRows() );
        }

        BOOST_CHECK_THROW (
        {
            setup.getGalerkin( numLevels );
        }, common::Exception );

        for ( IndexType i = 0; i < numLevels - 1; ++i )
        {
            const Matrix<ValueType>& galerkin0 = setup.getGalerkin( i );
            const Matrix<ValueType>& galerkin1 = setup.getGalerkin( i + 1 );
            const Matrix<ValueType>& interpolation = setup.getInterpolation( i  );
            const Matrix<ValueType>& restriction = setup.getRestriction( i  );

            // Interpolation: coarse ( galerkin1 ) -> fine ( galerkin0 )

            BOOST_CHECK_EQUAL( interpolation.getColDistribution(), galerkin1.getRowDistribution() );
            BOOST_CHECK_EQUAL( interpolation.getRowDistribution(), galerkin0.getRowDistribution() );

            BOOST_CHECK_EQUAL( restriction.getColDistribution(), interpolation.getRowDistribution() );
            BOOST_CHECK_EQUAL( restriction.getRowDistribution(), interpolation.getColDistribution() );
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();

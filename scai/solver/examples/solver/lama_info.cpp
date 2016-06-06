/**
 * @file lama_info.cpp
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
 * @brief Example program that prints internal information about LAMA factories
 * @author Thomas Brandes
 * @date 25.07.2015
 */

// Define levels for assertion, logging and tracing

#include "scai/lama.hpp"

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/solver/Solver.hpp>
#include <scai/solver/AMGSetup.hpp>

#include <scai/dmemo/NoCommunicator.hpp>

#include <scai/common/ContextType.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>

using namespace std;

void contextInfo()
{
    using namespace scai::hmemo;

    vector<Context::ContextType> values;  // supported context types

    Context::getCreateValues( values );

    cout << endl;
    cout << "Factory of Context: " << values.size() << " entries" << endl;
    cout << "=============================" << endl;
    cout << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "Registered values[" << i << "] = " << values[i] << endl;

        // ContextFactory has additional argument for device id

        ContextPtr context = Context::create( values[i], -1 );

        cout << "  Context: " << *context << endl;
        cout << "    Context->Memory: " << *context->getMemoryPtr() << endl;
        cout << "    Context->Host Memory: " << *context->getHostMemoryPtr() << endl;
    }

    cout << endl;
}

void communicatorInfo()
{
    using namespace scai::dmemo;

    vector<Communicator::CommunicatorKind> values;  // string is create type for the factory

    Communicator::getCreateValues( values );

    cout << endl;
    cout << "Factory of Communicator: " << values.size() << " entries" << endl;
    cout << "==================================" << endl;
    cout << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << values[i] << endl;
        CommunicatorPtr comm = Communicator::create( values[i]);
        cout << "    Communicator: " << *comm << endl;
    }

    cout << endl;
}

void matrixInfo()
{
    using namespace scai::lama;

    vector<MatrixCreateKeyType> keys;

    Matrix::getCreateValues( keys );

    cout << endl;
    cout << "Factory of Matrix: " << keys.size() << " entries" << endl;
    cout << "=============================" << endl;
    cout << endl;

    for ( size_t i = 0; i < keys.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << keys[i].first << ", " << keys[i].second << endl;

        scai::common::shared_ptr<Matrix> matrix ( Matrix::create( keys[i] ) );

        cout << "    Matrix: " << *matrix << endl;
    }

    cout << endl;
}

void vectorInfo()
{
    using namespace scai::lama;

    vector<VectorCreateKeyType> keys;

    Vector::getCreateValues( keys );

    cout << endl;
    cout << "Factory of Vector: " << keys.size() << " entries" << endl;
    cout << "=============================" << endl;
    cout << endl;

    for ( size_t i = 0; i < keys.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << keys[i].first << ", " << keys[i].second << endl;

        scai::common::shared_ptr<Vector> vector ( Vector::create( keys[i] ) );

        cout << "    Vector: " << *vector << endl;
    }

    cout << endl;
}

void solverInfo()
{
    using namespace scai::solver;

    vector<string> values;  // string is create type for the factory

    Solver::getCreateValues( values );

    cout << endl;
    cout << "Factory of Solver: " << values.size() << " entries" << endl;
    cout << "=============================" << endl;
    cout << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "   Registered values[" << i << "] = " << values[i] << endl;

        scai::common::shared_ptr<Solver> solver( Solver::create( values[i], "TestSolver" ) );

        cout << "      Solver: " << *solver << endl;
    }

    cout << endl;
}

void setupInfo()
{
    using namespace scai::lama;
    using namespace scai::solver;

    vector<string> values;  // string is create type for the factory

    AMGSetup::getCreateValues( values );

    cout << endl;
    cout << "Factory of AMG Setups: " << values.size() << " entries" << endl;
    cout << "================================" << endl;
    cout << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << values[i] << endl;

        scai::common::shared_ptr<AMGSetup> setup( AMGSetup::create( values[i] ) );

        cout << "    Setup: " << *setup << endl;
    }

    cout << endl;
}

void distributionInfo()
{
    using namespace scai::dmemo;

    vector<string> values;  // string is create type for the factory

    Distribution::getCreateValues( values );

    cout << endl;
    cout << "Factory of Distributions: " << values.size() << " entries" << endl;
    cout << "===================================" << endl;
    cout << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << values[i] << endl;

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // get the default one

        scai::common::shared_ptr<Distribution> dist( Distribution::getDistributionPtr( values[i], comm, 10, 1.0 ) );

        cout << "    Distribution: " << *dist << endl;
    }

    cout << endl;
}

int main( int /*argc */, char** /*argv*/ )
{
    std::string loadPath;

    if ( scai::common::Settings::getEnvironment( loadPath, "SCAI_LIBRARY_PATH" ) )
    {
        cout << "Load all module libraries in loadPath " << loadPath << endl;

        scai::common::LibModule::loadLibsByPath( loadPath.c_str() );
    }
    else
    {
        cout << "No additional modules loaded" << endl;
    }

    matrixInfo();
    vectorInfo();
    communicatorInfo();
    contextInfo();
    solverInfo();
    setupInfo();
    distributionInfo();
}


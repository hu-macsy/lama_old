/**
 * @file common/examples/DemoSettings.cpp
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
 * @brief Example of using the settings
 * @author Thomas Brandes
 * @date 27.01.2016
 */

#include <scai/common/Settings.hpp>

#include <iostream>

/* -----------------------------------------------------------------------------*/

using scai::common::Settings;

using namespace std;

int main( int argc, char const* argv[] )
{
    // print all SCAI values
    cout << "SCAI variables of environment:" << endl;
    cout << "==============================" << endl;
    Settings::printEnvironment( cout );
    cout << endl;

    Settings::parseArgs( argc, argv );

    if ( argc > 1 )
    {
        // there are some remaining arguments, shouldn't happpen
        cerr << "ERROR: only arguments like --SCAI_xxx are accepted." << endl;
        return -1;
    }

    cout << "SCAI variables of environment (after parsing command line args):" << endl;
    cout << "================================================================" << endl;
    Settings::printEnvironment( cout );
    cout << endl;

    string settingsFile;

    const char NAME[] = "scai";
    const int NODE_RANK = 3;

    if ( Settings::getEnvironment( settingsFile, "SCAI_SETTINGS" ) )
    {
        cout << "Read file " << settingsFile << " to find settings for node = " << NAME << ", rank = " << NODE_RANK << endl;
        Settings::readSettingsFile( settingsFile.c_str(), NAME, NODE_RANK );
        cout << "SCAI variables of environment (after reading settings):" << endl;
        cout << "================================================================" << endl;
        Settings::printEnvironment( cout );
        cout << endl;
    }

    // take the 4-th argument for comma separated value lists

    Settings::setRank( NODE_RANK );

    bool useMKL;

    if ( Settings::getEnvironment( useMKL, "SCAI_USE_MKL" ) )
    {
        cout << "useMKL = " << useMKL << endl;
    }
    else
    {
        cout << "SCAI_USE_MKL not defined" << endl;
    }

    int  device;

    if ( Settings::getEnvironment( device, "SCAI_DEVICE" ) )
    {
        cout << "device = " << device << endl;
    }
    else
    {
        cout << "SCAI_DEVICE not defined" << endl;
    }

    std::string solverName;

    if ( Settings::getEnvironment( solverName, "SCAI_SOLVER" ) )
    {
        cout << "solver = " << solverName << endl;
    }
    else
    {
        cout << "SCAI_SOLVER not defined" << endl;
    }
}


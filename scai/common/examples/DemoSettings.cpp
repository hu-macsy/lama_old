/**
 * @file common/examples/DemoSettings.cpp
 * @brief Example of using the settings
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
    Settings::printEnvironment();
    cout << endl;

    Settings::parseArgs( argc, argv );

    cout << "SCAI variables of environment (after parsing command line args):" << endl;
    cout << "================================================================" << endl;
    Settings::printEnvironment();
    cout << endl;

    // take the 4-th argument for comma separated value lists

    Settings::setRank( 3 );

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


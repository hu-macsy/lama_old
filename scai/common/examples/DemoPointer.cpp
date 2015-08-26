/**
 * @file common/examples/DemoPointer.cpp
 * @brief Example of pointer use 
 */

#include <scai/common/unique_ptr.hpp>

#include <iostream>

/* -----------------------------------------------------------------------------*/

int main()
{
    scai::common::unique_ptr<double> sum( new double );
    scai::common::scoped_array<double> vals ( new double[10] );

    *sum = 0.0;

    for ( int i = 0; i < 10; ++ i )
    {
       vals[i] = i;
    }

    for ( int i = 0; i < 10; ++ i )
    {
       *sum += vals[i];
    }

    std::cout << "Sum = " << *sum << std::endl;
}


/**
 * @file common/examples/TimePrecision.cpp
 *
 * @brief Check resolution of walltime
 *
 * @author Thomas Brandes
 * @date 10.12.2015
 */

#include <scai/common/Walltime.hpp>

#include <iostream>

using scai::common::Walltime;

int main()
{
   const int NITER = 4;

   for ( int i = 0; i < NITER; ++i )
   {
       long counter = 0;     // count changes of output value of Walltime.:get

       double t = Walltime::get();

       bool stop = false;   // set it to true after one second

       double tc = t;

       while ( !stop )
       {
          double t1 = Walltime::get();
    
          if ( t1 > tc )
          {
             // value has changed
    
             counter++;
             tc = t1;
          }
    
          stop = ( t1 - t ) >= 1.0;
       }
    
       std::cout << "Resolution: at least " << counter << " ticks per seconds" << std::endl;
    }
}

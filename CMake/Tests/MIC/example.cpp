/**
 * @file MIC/example.cpp
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Simple example program that might be used to test compilation for MIC device
 * @author Thomas Brandes
 * @date 12.06.2014
 */

#pragma offload_attribute (push,target(mic))
#include <iostream>
#include <omp.h>
#pragma offload_attribute (pop)

#include <cstdlib>
#include <string>
#include <stdio.h>

#include <offload.h> // to define _Offload_number_of_devices() and _Offload_get_device_number()

__attribute__ ((target(mic)))
  void hello_world_from_any();

__attribute__ ((target(mic)))
void testThreadCount();

__attribute__ ((target(mic)))
std::string get_device_name();

int main(int argc, char *argv[])
{
  hello_world_from_any();
  testThreadCount();

#pragma offload target(mic:0)
  {
    hello_world_from_any();
    testThreadCount();  
  }


#pragma offload target(mic:1)
  {
    hello_world_from_any();
    testThreadCount();  
  }

  exit(0);
}

std::string get_device_name()
{
 #ifdef __MIC__
  std::string device_name = "MIC ";
#else
  std::string device_name = "HOST";
#endif // __MIC__

 return device_name;
}

void hello_world_from_any()
{
  std::cout << "***** Hello by " << get_device_name() << " *****" << std::endl;
} 


void testThreadCount()
{
  int mic_rank = _Offload_get_device_number();
  int num_mics = _Offload_number_of_devices();

   #pragma omp parallel
   {
     int tid = omp_get_thread_num();
     #pragma omp critical
     {  
       std::cout << "Hello by thread id: " << tid << " --> " << get_device_name() << ": " << mic_rank << " / " << num_mics << std::endl;
     }
     #pragma omp barrier
     #pragma omp single
     {
       int thread_count = omp_get_num_threads();
       #pragma omp critical
       {
         std::cout << "Thread count: " << thread_count << std::endl;
         std::cout << "======================================" << std::endl;
       }
     }  
  }
}

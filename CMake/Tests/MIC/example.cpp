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

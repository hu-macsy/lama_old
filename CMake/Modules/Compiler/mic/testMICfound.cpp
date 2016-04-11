#pragma offload_attribute (push,target(mic))
#include <iostream>
#include <omp.h>
#pragma offload_attribute (pop)

#include <cstdlib>
#include <string>
#include <stdio.h>

#include <offload.h> // to define _Offload_number_of_devices() and _Offload_get_device_number()

__attribute__ ((target(mic)))
std::string get_device_name();

int main(int argc, char *argv[])
{
  std::cout << get_device_name() << std::endl;
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

}

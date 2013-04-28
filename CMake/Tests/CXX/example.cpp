// Simple example program that prints size of pointer (32-bit or 64-bit)

#include <iostream>

int main( int, char* [] )
{
    std::cout << "Size of (void*) = " << sizeof( void* ) << std::endl;

    return 0;
}

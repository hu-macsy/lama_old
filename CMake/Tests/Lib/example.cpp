// Simple example program that prints size of pointer (32-bit or 64-bit)

#include <iostream>
#include "Rectangle.hpp"

int main( int, char* [] )
{
    double x = 4.0;
    double y = 2.5;

    Rectangle r( x, y );

    std::cout << "Size of rectangle( x = " << x << ", y = " << y << " ) is " << r.getArea() << std::endl;

    return 0;
}

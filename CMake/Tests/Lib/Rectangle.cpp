
#define COMPILING_DLL

#include "Rectangle.hpp"

void Rectangle::setValues ( double _x, double _y )
{
    x = _x;
    y = _y;
}

Rectangle::Rectangle( double x, double y )
{
    setValues( x, y );
}

double Rectangle::getArea ( void )
{
    return x * y;
}

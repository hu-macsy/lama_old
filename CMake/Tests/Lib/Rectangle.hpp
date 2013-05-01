#include "config.hpp"

class DLL_IMPORTEXPORT Rectangle 
{
    double x;
    double y;

public:

    Rectangle( double, double );
    void setValues ( double, double );
    double getArea (void);
};

#include <cmath>


namespace scai {

namespace common {

float my_sqrt( const float& x )
{
	return sqrtf( x );
}

double my_sqrt( const double& x )
{
	return sqrt( x );
}

Complex<float> my_sqrt( const Complex<float>& x )
{
	return sqrt( x );
}

Complex<double> my_sqrt( const Complex<double>& x )
{
	return sqrt( x );
}

float my_fabs( const float& x )
{
	return fabsf( x );
}

double my_fabs( const double& x )
{
	return fabs( x );
}
Complex<float> my_fabs( const Complex<float>& x )
{
	return fabs( x );
}
Complex<double> my_fabs( const Complex<double>& x )
{
	return fabs( x );
}

} /* end namespace common */

} /* end namespace scai */

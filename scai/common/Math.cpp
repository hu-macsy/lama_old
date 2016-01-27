#include <scai/common/Math.hpp>

#include <cmath>

namespace scai {

namespace common {

float Math::sqrt( const float& x )
{
	return ::sqrtf( x );
}

double Math::sqrt( const double& x )
{
	return ::sqrt( x );
}

long double Math::sqrt( const long double& x )
{
	return ::sqrtl( x );
}

float Math::abs( const float& x )
{
	return ::fabsf( x );
}

double Math::abs( const double& x )
{
	return ::fabs( x );
}

long double Math::abs( const long double& x )
{
	return ::fabsl( x );
}

} /* end namespace common */

} /* end namespace scai */

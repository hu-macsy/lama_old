#include <scai/common/Math.hpp>

#include <cmath>

namespace scai {

namespace common {

MIC_CALLABLE_MEMBER float Math::sqrt( const float& x )
{
	return ::sqrtf( x );
}

MIC_CALLABLE_MEMBER double Math::sqrt( const double& x )
{
	return ::sqrt( x );
}

MIC_CALLABLE_MEMBER long double Math::sqrt( const long double& x )
{
	return ::sqrtl( x );
}

MIC_CALLABLE_MEMBER float Math::abs( const float& x )
{
	return ::fabsf( x );
}

MIC_CALLABLE_MEMBER double Math::abs( const double& x )
{
	return ::fabs( x );
}

MIC_CALLABLE_MEMBER long double Math::abs( const long double& x )
{
	return ::fabsl( x );
}

} /* end namespace common */

} /* end namespace scai */

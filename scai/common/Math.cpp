#include <scai/common/Math.hpp>

#include <cmath>
#include <cstdlib>

namespace scai {

namespace common {
// -------------------------------- sqrt -----------------------------
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

// -------------------------------- abs -----------------------------
MIC_CALLABLE_MEMBER int Math::abs( const int& x )
{
	return ::abs( x );
}

MIC_CALLABLE_MEMBER long Math::abs( const long& x )
{
	return ::labs( x );
}

MIC_CALLABLE_MEMBER long long Math::abs( const long long& x )
{
	return ::llabs( x );
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

// -------------------------------- conj -----------------------------
MIC_CALLABLE_MEMBER float Math::conj( const float& x )
{
	return x;
}

MIC_CALLABLE_MEMBER double Math::conj( const double& x )
{
	return x;
}

MIC_CALLABLE_MEMBER long double Math::conj( const long double& x )
{
	return x;
}

} /* end namespace common */

} /* end namespace scai */

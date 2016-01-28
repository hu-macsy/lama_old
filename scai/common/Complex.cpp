#include <scai/common/Complex.hpp>
#include <scai/common/Math.hpp>

#include <algorithm>

namespace scai {

namespace common {

template<>
long double Complex<long double>::metrikCuda( void ) const
{
	return Math::sqrt( real() * real() + imag() * imag() );
}

// ------------------ Math::sqrt --------------------------------
ComplexFloat Math::sqrt( const ComplexFloat& a )
{
    double x = a.real();
    double y = a.imag();

    if( x == 0.0f )
    {
    	double t = Math::sqrt( Math::abs( y ) / 2 );
        return ComplexFloat( t, y < 0.0f ? -t : t );
    }
    else
    {
    	double t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
    	double u = t / 2;
        return x > 0.0f ? ComplexFloat( u, y / t ) :
        		ComplexFloat( Math::abs( y ) / t, y < 0.0f ? -u : u );
    }
}

ComplexDouble Math::sqrt( const ComplexDouble& a )
{
    double x = a.real();
    double y = a.imag();

    if( x == 0.0 )
    {
    	double t = Math::sqrt( Math::abs( y ) / 2 );
        return ComplexDouble( t, y < 0.0 ? -t : t );
    }
    else
    {
    	double t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
    	double u = t / 2;
        return x > 0.0 ? ComplexDouble( u, y / t ) :
        		ComplexDouble( Math::abs( y ) / t, y < 0.0 ? -u : u );
    }
}

ComplexLongDouble Math::sqrt( const ComplexLongDouble& a )
{
	long double x = a.real();
	long double y = a.imag();

    if( x == 0.0l )
    {
    	long double t = Math::sqrt( Math::abs( y ) / 2 );
        return ComplexLongDouble( t, y < 0.0l ? -t : t );
    }
    else
    {
    	long double t = Math::sqrt( 2 * ( abs( a ) + Math::abs( x ) ) );
    	long double u = t / 2;
        return x > 0.0l ? ComplexLongDouble( u, y / t ) :
        		ComplexLongDouble( Math::abs( y ) / t, y < 0.0l ? -u : u );
    }
}

// ------------------ Math::abs --------------------------------
float Math::abs( const ComplexFloat& a )
{
	float x = a.real();
	float y = a.imag();
	const float s = std::max( Math::abs(x), Math::abs(y));
	if (s == 0.0f)
	{
		return s;
	}
	x /= s;
	y /= s;
	return s * Math::sqrt( x * x + y * y);
}

double Math::abs( const ComplexDouble& a )
{
	double x = a.real();
	double y = a.imag();
	const double s = std::max( Math::abs(x), Math::abs(y) );
	if (s == 0.0)
	{
		return s;
	}
	x /= s;
	y /= s;
	return s * Math::sqrt( x * x + y * y);
}

long double Math::abs( const ComplexLongDouble& a )
{
	long double x = a.real();
	long double y = a.imag();
	const long double s = std::max( Math::abs(x), Math::abs(y) );
	if (s == 0.0l)
	{
		return s;
	}
	x /= s;
	y /= s;
	return s * Math::sqrt( x * x + y * y);
}

// ------------------ Math::conj --------------------------------
ComplexFloat Math::conj( const ComplexFloat& a )
{
	return ComplexFloat( a.real(), -a.imag() );
}
ComplexDouble Math::conj( const ComplexDouble& a )
{
	return ComplexDouble( a.real(), -a.imag() );
}
ComplexLongDouble Math::conj( const ComplexLongDouble& a )
{
	return ComplexLongDouble( a.real(), -a.imag() );
}

} /* end namespace common */

} /* end namespace scai */

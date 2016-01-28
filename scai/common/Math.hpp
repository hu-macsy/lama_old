/*
 * math.hpp
 *
 *  Created on: Jan 26, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/mic/MICCallable.hpp>
#include <scai/common/cuda/CUDACallable.hpp>

namespace scai {

namespace common {


#ifdef SCAI_COMPLEX_SUPPORTED
	template<typename ValueType> class Complex;
#endif

struct Math
{
	/** Square root function for ValueType
	 *
	 *  In contrary to the routine of cmath it will be possible to
	 *  use always the same name for the routine.
	 *
	 *  \code
	 *    ValueType x = sqrt ( y );                          // might not work always correctly
	 *    ValueType x = TypeTraits<ValueType>::sqrt ( y );   // this is guaranteed to work
	 *  \endcode
	 */
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sqrt( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sqrt( const double& x );

	static MIC_CALLABLE_MEMBER long double sqrt( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sqrt( const Complex<float>& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sqrt( const Complex<double>& x );

	static MIC_CALLABLE_MEMBER Complex<long double> sqrt( const Complex<long double>& x );
#endif

    /** Absolute value function for ValueType
     *
     *  In contrary to the routine of cmath it will be possible to
     *  use always the same name for the routine.
     */
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER int abs( const int& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long abs( const long& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long long abs( const long long& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const double& x );

	static MIC_CALLABLE_MEMBER long double abs( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const Complex<float>& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const Complex<double>& x );

	static MIC_CALLABLE_MEMBER long double abs( const Complex<long double>& x );
#endif

	/*
	 * Computes the conjugated value of a given value
	 */
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float conj( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double conj( const double& x );

	static MIC_CALLABLE_MEMBER long double conj( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> conj( const Complex<float>& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> conj( const Complex<double>& x );

	static MIC_CALLABLE_MEMBER Complex<long double> conj( const Complex<long double>& x );
#endif

};

} /* end namespace common */

} /* end namespace scai */

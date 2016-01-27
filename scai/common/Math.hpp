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
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sqrt( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sqrt( const double& x );

	static long double sqrt( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<float> sqrt( const Complex<float>& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER Complex<double> sqrt( const Complex<double>& x );

	static Complex<long double> sqrt( const Complex<long double>& x );
#endif

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const double& x );

	static long double abs( const long double& x );

#ifdef SCAI_COMPLEX_SUPPORTED
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float abs( const Complex<float>& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double abs( const Complex<double>& x );

	static long double abs( const Complex<long double>& x );
#endif
};

} /* end namespace common */

} /* end namespace scai */

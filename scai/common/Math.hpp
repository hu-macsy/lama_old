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

struct Math
{
	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float sqrt( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double sqrt( const double& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long double sqrt( const long double& x );


	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER float fabs( const float& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER double fabs( const double& x );

	static MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER long double fabs( const long double& x );
};

} /* end namespace common */

} /* end namespace scai */

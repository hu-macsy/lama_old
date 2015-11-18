#pragma once

#include <mkl.h>

class MKLUtils
{
public:
#ifdef __INTEL_OFFLOAD
	__declspec( target(mic) )
#endif
	static inline MKL_Complex8* cast(ComplexFloat* in)
	{
		return reinterpret_cast<MKL_Complex8*>( in );
	}

#ifdef __INTEL_OFFLOAD
	__declspec( target(mic) )
#endif
	static inline const MKL_Complex8* cast(const ComplexFloat* in)
	{
		return reinterpret_cast<const MKL_Complex8*>( in );
	}

#ifdef __INTEL_OFFLOAD
	__declspec( target(mic) )
#endif
	static inline MKL_Complex16* cast(ComplexDouble* in)
	{
		return reinterpret_cast<MKL_Complex16*>( in );
	}

#ifdef __INTEL_OFFLOAD
	__declspec( target(mic) )
#endif
	static inline const MKL_Complex16* cast(const ComplexDouble* in)
	{
		return reinterpret_cast<const MKL_Complex16*>( in );
	}
};

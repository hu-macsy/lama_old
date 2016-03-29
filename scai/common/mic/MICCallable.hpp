
#pragma once

//#if defined( __INTEL_OFFLOAD ) && defined( __MIC__ )
#if defined( __INTEL_OFFLOAD )
    #define MIC_CALLABLE_MEMBER __declspec( target(mic) )
#else
    #define MIC_CALLABLE_MEMBER
#endif

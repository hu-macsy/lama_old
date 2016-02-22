
#pragma once

#ifdef __INTEL_OFFLOAD
    #define MIC_CALLABLE_MEMBER __declspec( target(mic) )
#else
    #define MIC_CALLABLE_MEMBER
#endif

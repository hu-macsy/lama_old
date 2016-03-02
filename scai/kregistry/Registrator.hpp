#pragma once

#include <scai/kregistry/KernelRegistry.hpp>

#define SCAI_DECLARE_REGISTRATOR( name, ... )                                                            \
    template< __VA_ARGS__ >                                                            \
    struct name                                                                      \
    {                                                                                       \
        static void initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag );       \
    };


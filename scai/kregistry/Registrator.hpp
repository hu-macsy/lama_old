#pragma once

#include <scai/kregistry/KernelRegistry.hpp>

#define SCAI_DECLARE_REGISTRATOR( ValueType )                                                            \
    template<typename ValueType>                                                            \
    struct Registrator                                                                      \
    {                                                                                       \
        static void initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag );       \
    };


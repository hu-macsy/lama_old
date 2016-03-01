#pragma once

#define SCAI_DECLARE_REGISTRATOR( ValueType )                                                            \
    template<typename ValueType>                                                            \
    struct Registrator                                                                      \
    {                                                                                       \
        static void initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag );       \
    };


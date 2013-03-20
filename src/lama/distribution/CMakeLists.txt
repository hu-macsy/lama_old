### set include directories ###

include_directories ( . )

set ( CXXSOURCES
         ${CXXSOURCES}
         distribution/BlockDistribution.cpp
         distribution/BlockDistributionManager.cpp
         distribution/CyclicDistribution.cpp
         distribution/CyclicDistributionManager.cpp
         distribution/Distribution.cpp
         distribution/GenBlockDistribution.cpp
         distribution/GeneralDistribution.cpp
         distribution/NoDistribution.cpp
         distribution/NoDistributionManager.cpp
         distribution/Halo.cpp
         distribution/HaloBuilder.cpp
         distribution/Redistributor.cpp
         
         PARENT_SCOPE
    )

### set include directories ### 

include_directories ( . )


### set source files ###

set ( CXXSOURCES
        ${CXXSOURCES}
        norm/L1Norm.cpp
        norm/L2Norm.cpp
        norm/MaxNorm.cpp
        norm/Norm.cpp
        PARENT_SCOPE
    )

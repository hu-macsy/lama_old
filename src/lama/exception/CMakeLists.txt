### set include directories ### 

include_directories ( . )

### set source files ###

set ( CXXSOURCES
        ${CXXSOURCES}
        exception/Exception.cpp 
        # add new exceptions cxx files here
        PARENT_SCOPE
    )

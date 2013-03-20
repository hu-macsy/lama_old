### set include directories ### 

include_directories ( . )


### set source files ###

set ( CXXSOURCES
        ${CXXSOURCES}
        storage/COOStorage.cpp
        storage/CSRStorage.cpp
        storage/DenseStorage.cpp
        storage/DIAStorage.cpp
        storage/ELLStorage.cpp
        storage/JDSStorage.cpp
        storage/MatrixStorage.cpp
        storage/StorageMethods.cpp
        storage/SparseAssemblyStorage.cpp
        PARENT_SCOPE
    )

### set include directories ### 

include_directories ( . )


### set source files ###

set ( CXXSOURCES
        ${CXXSOURCES}
      	matrix/COOSparseMatrix.cpp
      	matrix/CSRSparseMatrix.cpp
      	matrix/DenseMatrix.cpp
      	matrix/DIASparseMatrix.cpp
      	matrix/ELLSparseMatrix.cpp
      	matrix/JDSSparseMatrix.cpp
      	matrix/Matrix.cpp
      	matrix/SparseMatrix.cpp
      	#matrix/MetaMatrix.cpp
      	#matrix/SimpleStorageStrategy.cpp
        PARENT_SCOPE
    )

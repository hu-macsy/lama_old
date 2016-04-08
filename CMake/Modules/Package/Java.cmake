if    ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
	# Enable Java Compilation
	# this will set at least CMAKE_Java_COMPILER CMAKE_Java_ARCHIVE 
    find_package( Java )
endif ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
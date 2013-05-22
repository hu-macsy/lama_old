###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief logging/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
 # $Id$
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        GenLogger
        GenLoggerCreator
        Level
        Logger
        LoggerProvider
        SourceLocation
    )

set ( HEADERS
        AbstractLoggerCreator
        config
        logging
        NonCopyable
    )
    
lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )

### add definitions ###
# GenLogger uses also logging
add_definitions ( -DLAMA_LOG_LEVEL_DEBUG )

if ( WIN32 )
	add_definitions ( -DLOG4LAMA_COMPILING_DLL )
endif ( WIN32 )

### add library and put path to search path ###
add_library ( log4lama SHARED ${CXX_SOURCES} )
include_directories ( . )

### install ###
install ( TARGETS log4lama DESTINATION lib )

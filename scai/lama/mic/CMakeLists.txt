###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief lama/mic/CMakeLists.txt
 # @author Thomas Brandes
 # @date 01.07.2013
 # @since 1.0.1
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        MICContext
        MICSyncToken
        MICContextManager
        MICUtils
        MICCSRUtils
        MICMKLCSRUtils
        MICELLUtils
        MICJDSUtils
        MICDIAUtils
        MICCOOUtils
        MICBLAS1
        MICBLAS2
        MICBLAS3
    )

lama_classes ( ${CLASSES} )

add_library( amamic ${LAMA_LIBRARY_TYPE} ${CXX_SOURCES} )

### target link libraries ###


set ( MKLROOT /cluster/intel/composer_xe_2013.1.117/mkl )

set_target_properties( amamic PROPERTIES LINK_FLAGS "-offload-option,mic,compiler,\"-L${MKLROOT}/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group\"" )

### install library
install ( TARGETS amamic DESTINATION lib )

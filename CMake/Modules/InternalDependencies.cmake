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
 # @brief Central defenition of internal dependencies between sub projects
 # @author Lauretta Schubert
 # @date 17.08.2015
 # @since 2.0.0
###

set ( SCAI_COMMON_INTERNAL_DEPS )
set ( SCAI_LOGGING_INTERNAL_DEPS scai_common )
set ( SCAI_TRACING_INTERNAL_DEPS scai_common scai_logging )
set ( SCAI_TASKING_INTERNAL_DEPS scai_common scai_logging scai_tracing )
#set ( SCAI_KERNEL_INTERNAL_DEPS scai_common scai_logging scai_tracing )
set ( SCAI_MEMORY_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking )
set ( SCAI_LAMA_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking scai_memory ) #scai_kernel )

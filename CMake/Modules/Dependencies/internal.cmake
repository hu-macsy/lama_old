###
 # @file dependendies/internal.cmake
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

## CAUTION: define internal dependencies considering their dependencies !!!
## use reverse linking order (for static linking)

set ( SCAI_COMMON_INTERNAL_DEPS ) # 1
set ( SCAI_LOGGING_INTERNAL_DEPS      scai_common ) # 2
set ( SCAI_TRACING_INTERNAL_DEPS      scai_common scai_logging ) # 3
set ( SCAI_TASKING_INTERNAL_DEPS      scai_common scai_logging scai_tracing ) # 4
set ( SCAI_HMEMO_INTERNAL_DEPS        scai_common scai_logging scai_tracing scai_tasking ) # 5
set ( SCAI_KREGISTRY_INTERNAL_DEPS    scai_common scai_logging ) # 6
set ( SCAI_BLASKERNEL_INTERNAL_DEPS   scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry ) # 7
set ( SCAI_UTILSKERNEL_INTERNAL_DEPS  scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry ) # 8
set ( SCAI_SPARSEKERNEL_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_utilskernel ) # 9
set ( SCAI_DMEMO_INTERNAL_DEPS        scai_common scai_logging scai_tracing scai_tasking scai_hmemo ) # 10
set ( SCAI_LAMA_INTERNAL_DEPS         scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo ) # 11
set ( SCAI_SOLVER_INTERNAL_DEPS       scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo scai_lama ) #12

## head project containing all sub projects
set ( LAMA_ALL_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo scai_lama scai_solver )

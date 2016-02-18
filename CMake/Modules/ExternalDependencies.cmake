###
 # @file ExternalDependencies.cmake
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
 # @brief Central defenition of external dependencies between sub projects
 # @author Lauretta Schubert
 # @date 17.08.2015
 # @since 2.0.0
###

## attention OpenMP should be before SCAI_BLAS !!!
## need to get rid of boost for boost_repeat in common

set ( SCAI_COMMON_EXTERNAL_DEPS dl Thread Boost OpenMP CUDA ) # 1
set ( SCAI_LOGGING_EXTERNAL_DEPS ) # 2
set ( SCAI_TRACING_EXTERNAL_DEPS ) # 3
set ( SCAI_TASKING_EXTERNAL_DEPS ) # 4
set ( SCAI_HMEMO_EXTERNAL_DEPS                         CUDA MIC ) # 5
set ( SCAI_KREGISTRY_EXTERNAL_DEPS ) # 6
set ( SCAI_BLASKERNEL_EXTERNAL_DEPS             OpenMP CUDA MIC SCAI_BLAS ) # 7
set ( SCAI_DMEMO_EXTERNAL_DEPS                  OpenMP                    MPI GPI GraphPartitioning ) # 8
set ( SCAI_LAMA_EXTERNAL_DEPS                   OpenMP CUDA MIC SCAI_BLAS ) # 9
set ( SCAI_SOLVER_EXTERNAL_DEPS                                 SCAI_BLAS ) # 10

set ( SCAI_EXTERNAL_DEPS ${SCAI_COMMON_EXTERNAL_DEPS} ${SCAI_LOGGING_EXTERNAL_DEPS} ${SCAI_TRACING_EXTERNAL_DEPS}
						 ${SCAI_TASKING_EXTERNAL_DEPS} ${SCAI_HMEMO_EXTERNAL_DEPS} ${SCAI_KREGISTRY_EXTERNAL_DEPS}
						 ${SCAI_DMEMO_EXTERNAL_DEPS} ${SCAI_BLASKERNEL_EXTERNAL_DEPS} ${SCAI_LAMA_EXTERNAL_DEPS}
						 ${SCAI_SOLVER_EXTERNAL_DEPS} )

list ( REMOVE_DUPLICATES SCAI_EXTERNAL_DEPS )

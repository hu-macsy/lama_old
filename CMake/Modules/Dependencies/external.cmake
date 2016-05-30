###
 # @file Dependencies/external.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Central defenition of external dependencies between sub projects
 # @author Lauretta Schubert
 # @date 17.08.2015
###

## attention OpenMP should be before SCAI_BLAS !!!
## need to get rid of boost for boost_repeat in common

#                                     dl Java OpenMP CUDA MIC  SCAI_BLAS MPI GPI GraphPartitioning
#                                    -------------------------------------------------------------
set ( SCAI_COMMON_EXTERNAL_DEPS       dl      OpenMP CUDA MIC                                      ) # 1
set ( SCAI_LOGGING_EXTERNAL_DEPS                                                                   ) # 2
set ( SCAI_TRACING_EXTERNAL_DEPS         Java                                                      ) # 3
set ( SCAI_TASKING_EXTERNAL_DEPS              CUDA MIC                                             ) # 4
set ( SCAI_HMEMO_EXTERNAL_DEPS                OpenMP CUDA MIC                                      ) # 5
set ( SCAI_KREGISTRY_EXTERNAL_DEPS                                                                 ) # 6
set ( SCAI_BLASKERNEL_EXTERNAL_DEPS           OpenMP CUDA MIC SCAI_BLAS                            ) # 7
set ( SCAI_UTILSKERNEL_EXTERNAL_DEPS          OpenMP CUDA MIC                                      ) # 8
set ( SCAI_SPARSEKERNEL_EXTERNAL_DEPS         OpenMP CUDA MIC SCAI_BLAS                            ) # 9
set ( SCAI_DMEMO_EXTERNAL_DEPS                OpenMP                     MPI GPI GraphPartitioning ) # 10
set ( SCAI_LAMA_EXTERNAL_DEPS                                                                      ) # 11
set ( SCAI_SOLVER_EXTERNAL_DEPS                                                                    ) # 12

set ( LAMA_ALL_EXTERNAL_DEPS ${SCAI_COMMON_EXTERNAL_DEPS} ${SCAI_LOGGING_EXTERNAL_DEPS} ${SCAI_TRACING_EXTERNAL_DEPS}
						 	 ${SCAI_TASKING_EXTERNAL_DEPS} ${SCAI_HMEMO_EXTERNAL_DEPS} ${SCAI_KREGISTRY_EXTERNAL_DEPS}
						 	 ${SCAI_BLASKERNEL_EXTERNAL_DEPS} ${SCAI_UTILSKERNEL_EXTERNAL_DEPS}
						 	 ${SCAI_SPARSEKERNEL_EXTERNAL_DEPS} ${SCAI_DMEMO_EXTERNAL_DEPS} ${SCAI_LAMA_EXTERNAL_DEPS}
						 	 ${SCAI_SOLVER_EXTERNAL_DEPS} )

list ( REMOVE_DUPLICATES LAMA_ALL_EXTERNAL_DEPS )

###
 # @file Dependencies/internal.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 # @brief Central defenition of internal dependencies between sub projects
 # @author Lauretta Schubert
 # @date 17.08.2015
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
set ( SCAI_UTILSKERNEL_INTERNAL_DEPS  scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel ) # 8
set ( SCAI_SPARSEKERNEL_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_utilskernel ) # 9
set ( SCAI_DMEMO_INTERNAL_DEPS        scai_common scai_logging scai_tracing scai_tasking scai_hmemo ) # 10
set ( SCAI_LAMA_INTERNAL_DEPS         scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo ) # 11
set ( SCAI_SOLVER_INTERNAL_DEPS       scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo scai_lama ) #12

## head project containing all sub projects
set ( LAMA_ALL_INTERNAL_DEPS scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo scai_lama scai_solver )

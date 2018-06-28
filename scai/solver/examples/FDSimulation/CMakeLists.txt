###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration file for scai_lama examples
 # @author Vanessa Wolff
 # @date 24.05.2016
###

scai_add_example( EXECUTABLE FDSimulation.exe
                  FILES      FDSimulation.cpp )

# seismogram.gnuplot: copy only to build directory
configure_file( seismogram.gnuplot seismogram.gnuplot COPYONLY )
# Configuration.txt: copy only to build directory
configure_file( Configuration.txt Configuration.txt COPYONLY )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# generate build file in example directory
scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

## install examples sources and additional files
install ( FILES       ${EXAMPLE_FILES} seismogram.gnuplot Configuration.txt Configuration.hpp
          DESTINATION ${INSTALL_EXAMPLE_DIR} )

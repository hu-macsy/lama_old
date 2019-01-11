###
 # @file FindSphinx.cmake
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
 # @brief Find Sphinx
 # @author Jan Ecker
 # @date 24.07.2015
###

set (_Sphinx_REQUIRED_VARS)

# ----------------------------------------------------------------------------
# initialize search
if (NOT Sphinx_DIR)
  if (NOT $ENV{Sphinx_DIR} STREQUAL "")
    set (Sphinx_DIR "$ENV{Sphinx_DIR}" CACHE PATH "Installation prefix of Sphinx (docutils)." FORCE)
  else ()
    set (Sphinx_DIR "$ENV{SPHINX_DIR}" CACHE PATH "Installation prefix of Sphinx (docutils)." FORCE)
  endif ()
endif ()

# ----------------------------------------------------------------------------
# default components to look for
if (NOT Sphinx_FIND_COMPONENTS)
  set (Sphinx_FIND_COMPONENTS "build" "apidoc")
elseif (NOT Sphinx_FIND_COMPONENTS MATCHES "^(build|apidoc)$")
  message (FATAL_ERROR "Invalid Sphinx component in: ${Sphinx_FIND_COMPONENTS}")
endif ()

# ----------------------------------------------------------------------------
# find components, i.e., build tools
foreach (_Sphinx_TOOL IN LISTS Sphinx_FIND_COMPONENTS)
  if (Sphinx_DIR)
    find_program (
      Sphinx-${_Sphinx_TOOL}_EXECUTABLE
       NAMES         sphinx-${_Sphinx_TOOL} sphinx-${_Sphinx_TOOL}.py
      HINTS         "${Sphinx_DIR}"
      PATH_SUFFIXES bin
      DOC           "The sphinx-${_Sphinx_TOOL} Python script."
      NO_DEFAULT_PATH
    )
  else ()
    find_program (
      Sphinx-${_Sphinx_TOOL}_EXECUTABLE
      NAMES sphinx-${_Sphinx_TOOL} sphinx-${_Sphinx_TOOL}.py
      DOC   "The sphinx-${_Sphinx_TOOL} Python script."
    )
  endif ()
  mark_as_advanced (Sphinx-${_Sphinx_TOOL}_EXECUTABLE)
  list (APPEND _Sphinx_REQUIRED_VARS Sphinx-${_Sphinx_TOOL}_EXECUTABLE)
endforeach ()

# ----------------------------------------------------------------------------
# determine Python executable used by Sphinx
if (Sphinx-build_EXECUTABLE)
  # extract python executable from shebang of sphinx-build
  find_package (PythonInterp QUIET)
  set (Sphinx_PYTHON_EXECUTABLE "${PYTHON_EXECUTABLE}")
  set (Sphinx_PYTHON_OPTIONS)
  file (STRINGS "${Sphinx-build_EXECUTABLE}" FIRST_LINE LIMIT_COUNT 1)
  if (FIRST_LINE MATCHES "^#!(.*/python.*)") # does not match "#!/usr/bin/env python" !
    string (REGEX REPLACE "^ +| +$" "" Sphinx_PYTHON_EXECUTABLE "${CMAKE_MATCH_1}")
    if (Sphinx_PYTHON_EXECUTABLE MATCHES "([^ ]+) (.*)")
      set (Sphinx_PYTHON_EXECUTABLE "${CMAKE_MATCH_1}")
      string (REGEX REPLACE " +" ";" Sphinx_PYTHON_OPTIONS "${CMAKE_MATCH_2}")
    endif ()
  endif ()
  # this is done to avoid problems with multiple Python versions being installed
  # remember: CMake command if(STR EQUAL STR) is bad and may cause many troubles !
  string (REGEX REPLACE "([.+*?^$])" "\\\\\\1" _Sphinx_PYTHON_EXECUTABLE_RE "${PYTHON_EXECUTABLE}")
  list (FIND Sphinx_PYTHON_OPTIONS -E IDX)
  if (IDX EQUAL -1 AND NOT Sphinx_PYTHON_EXECUTABLE MATCHES "^${_Sphinx_PYTHON_EXECUTABLE_RE}$")
    list (INSERT Sphinx_PYTHON_OPTIONS 0 -E)
  endif ()
  unset (_Sphinx_PYTHON_EXECUTABLE_RE)
endif ()

# ----------------------------------------------------------------------------
# determine Sphinx version
if (Sphinx-build_EXECUTABLE)
  # intentionally use invalid -h option here as the help that is shown then
  # will include the Sphinx version information
  if (Sphinx_PYTHON_EXECUTABLE)
    execute_process (
      COMMAND "${Sphinx_PYTHON_EXECUTABLE}" ${Sphinx_PYTHON_OPTIONS} "${Sphinx-build_EXECUTABLE}" -h
      OUTPUT_VARIABLE _Sphinx_VERSION
      ERROR_VARIABLE  _Sphinx_VERSION
    )
  elseif (UNIX)
    execute_process (
      COMMAND "${Sphinx-build_EXECUTABLE}" -h
      OUTPUT_VARIABLE _Sphinx_VERSION
      ERROR_VARIABLE  _Sphinx_VERSION
    )
  endif ()
  if (_Sphinx_VERSION MATCHES "Sphinx v([0-9]+\\.[0-9]+\\.[0-9]+)")
    set (Sphinx_VERSION_STRING "${CMAKE_MATCH_1}")
    string (REPLACE "." ";" _Sphinx_VERSION "${Sphinx_VERSION_STRING}")
    list(GET _Sphinx_VERSION 0 Sphinx_VERSION_MAJOR)
    list(GET _Sphinx_VERSION 1 Sphinx_VERSION_MINOR)
    list(GET _Sphinx_VERSION 2 Sphinx_VERSION_PATCH)
    if (Sphinx_VERSION_PATCH EQUAL 0)
      string (REGEX REPLACE "\\.0$" "" Sphinx_VERSION_STRING "${Sphinx_VERSION_STRING}")
    endif ()
  endif()
endif ()

# ----------------------------------------------------------------------------
# compatibility with FindPythonInterp.cmake and FindPerl.cmake
set (SPHINX_EXECUTABLE "${Sphinx-build_EXECUTABLE}")

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set SPHINX_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (
  Sphinx
  REQUIRED_VARS
    ${_Sphinx_REQUIRED_VARS}
  VERSION_VAR
    Sphinx_VERSION_STRING
)

# ----------------------------------------------------------------------------
# set Sphinx_DIR
if (NOT Sphinx_DIR AND Sphinx-build_EXECUTABLE)
  get_filename_component (Sphinx_DIR "${Sphinx-build_EXECUTABLE}" PATH)
  string (REGEX REPLACE "/bin/?" "" Sphinx_DIR "${Sphinx_DIR}")
  set (Sphinx_DIR "${Sphinx_DIR}" CACHE PATH "Installation directory of Sphinx tools." FORCE)
endif ()

unset (_Sphinx_VERSION)
unset (_Sphinx_REQUIRED_VARS)
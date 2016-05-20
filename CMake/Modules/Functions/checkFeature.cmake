# Taken from http://ghulbus-inc.de/projects/stuff/cmake_cxx11.zip
# linked by http://pageant.ghulbus.eu/?p=664
#
# Original script by Rolf Eike Beer
# Modifications by Andreas Weis
#
# Adapted for LAMA by Lauretta Schubert
#
# for c++11 test example files see: CMake/Modules/Compiler/c++11/c++11-test-*.cpp

MACRO ( checkFeature FEATURE_NAME FEATURE_NUMBER RESULT_VAR COMPILE_FLAG )

	IF ( NOT DEFINED ${RESULT_VAR} )
		SET( _bindir "${CMAKE_CURRENT_BINARY_DIR}/c++11/cxx11/cxx11_${FEATURE_NAME}" )

		IF ( ${FEATURE_NUMBER} )
			SET ( _SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11/c++11-test-${FEATURE_NAME}-N${FEATURE_NUMBER} )
			SET ( _LOG_NAME "\"${FEATURE_NAME}\" (N${FEATURE_NUMBER})" )
		ELSE ( ${FEATURE_NUMBER} )
			SET ( _SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11/c++11-test-${FEATURE_NAME} )
			SET ( _LOG_NAME "\"${FEATURE_NAME}\"" )
		ENDIF ( ${FEATURE_NUMBER} )
		MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME}" )

		SET ( _SRCFILE "${_SRCFILE_BASE}.cpp" )
		SET ( _SRCFILE_FAIL "${_SRCFILE_BASE}_fail.cpp" )
		SET ( _SRCFILE_FAIL_COMPILE "${_SRCFILE_BASE}_fail_compile.cpp" )

		if    ( "${COMPILE_FLAG}" STREQUAL "" )
			set ( COMPILE_FLAG_OPTION "" )
		else  ( "${COMPILE_FLAG}" STREQUAL "" )
			#set ( COMPILE_FLAG_OPTION "CMAKE_FLAGS -DCMAKE_CXX_FLAGS=${COMPILE_FLAG}" )
			#set ( COMPILE_FLAG_OPTION "COMPILE_DEFINITIONS ${COMPILE_FLAG}" )
		endif ( "${COMPILE_FLAG}" STREQUAL "" )

		message ( STATUS "COMPILE_FLAG_OPTION: ${COMPILE_FLAG_OPTION}" )

		IF ( CROSS_COMPILING )
			try_compile ( ${RESULT_VAR} "${_bindir}" "${_SRCFILE}" )
			IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
				try_compile ( ${RESULT_VAR} "${_bindir}_fail" "${_SRCFILE_FAIL}" )
			ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
		ELSE ( CROSS_COMPILING )
			try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}" "${_SRCFILE}" ${COMPILE_FLAG_OPTION} )
			#try_compile ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}" "${_SRCFILE}" ${COMPILE_FLAG_OPTION} )

			#message ( STATUS "COMPILE_OUTPUT_VARIABLE ${COMP_VAR} " )
			#message ( STATUS )
			#message ( STATUS "RUN_OUTPUT_VAR ${RUN_VAR} " )


			IF ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
				SET ( ${RESULT_VAR} TRUE )
			ELSE ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
				SET ( ${RESULT_VAR} FALSE )
			ENDIF ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
			IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
				try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}_fail" "${_SRCFILE_FAIL}" )
				IF ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
					SET ( ${RESULT_VAR} TRUE )
				ELSE ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
					SET ( ${RESULT_VAR} FALSE )
				ENDIF ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
			ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
		ENDIF ( CROSS_COMPILING )

		IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE} )
			try_compile ( _TMP_RESULT "${_bindir}_fail_compile" "${_SRCFILE_FAIL_COMPILE}" )
			IF ( _TMP_RESULT )
				SET ( ${RESULT_VAR} FALSE )
			ELSE ( _TMP_RESULT )
				SET ( ${RESULT_VAR} TRUE )
			ENDIF ( _TMP_RESULT )
		ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE} )

		IF ( ${RESULT_VAR} )
			MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME} -- works" )
			#LIST ( APPEND CXX11_FEATURE_LIST ${RESULT_VAR} )
		ELSE ( ${RESULT_VAR} )
			MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME} -- not supported" )
			message ( STATUS "compile ${_COMPILE_RESULT_VAR}" )
			message ( STATUS "run ${_RUN_RESULT_VAR}" )
		ENDIF ( ${RESULT_VAR} )

		SET ( ${RESULT_VAR} ${${RESULT_VAR}} CACHE INTERNAL "C++11 support for ${_LOG_NAME}" )
	ENDIF ( NOT DEFINED ${RESULT_VAR} )

ENDMACRO ( checkFeature FEATURE_NAME FEATURE_NUMBER RESULT_VAR )

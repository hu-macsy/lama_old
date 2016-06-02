# Taken from http://ghulbus-inc.de/projects/stuff/cmake_cxx11.zip
# linked by http://pageant.ghulbus.eu/?p=664
#
# Original script by Rolf Eike Beer
# Modifications by Andreas Weis
#
# Adapted for LAMA by Lauretta Schubert
#
# for c++11 test example files see: CMake/Modules/Compiler/c++11/c++11-test-*.cpp


# COMPILE_FLAG is ignored, does not work with try_compile/try_run
# do the following instead:
# save old DCMAKE_CXX_FLAGS, set them to new value and restore old afterwards
MACRO ( checkFeature FEATURE_NAME FEATURE_NUMBER RESULT_VAR COMPILE_FLAG LINK_LIB )

	IF ( NOT DEFINED ${RESULT_VAR} )
		SET( _bindir "${CMAKE_CURRENT_BINARY_DIR}/c++11/cxx11/cxx11_${FEATURE_NAME}" )

		IF ( ${FEATURE_NUMBER} )
			SET ( _SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11/c++11-test-${FEATURE_NAME}-N${FEATURE_NUMBER} )
			SET ( _LOG_NAME "\"${FEATURE_NAME}\" (N${FEATURE_NUMBER})" )
		ELSE ( ${FEATURE_NUMBER} )
			SET ( _SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11/c++11-test-${FEATURE_NAME} )
			SET ( _LOG_NAME "\"${FEATURE_NAME}\"" )
		ENDIF ( ${FEATURE_NUMBER} )
		#MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME}" )

		SET ( _SRCFILE "${_SRCFILE_BASE}.cpp" )
		SET ( _SRCFILE_FAIL "${_SRCFILE_BASE}_fail.cpp" )
		SET ( _SRCFILE_FAIL_COMPILE "${_SRCFILE_BASE}_fail_compile.cpp" )

		set ( OPTION_ON FALSE )

		#if    ( "${COMPILE_FLAG}" STREQUAL "" )
			set ( COMPILE_FLAG_OPTION "" )
		#else  ( "${COMPILE_FLAG}" STREQUAL "" )
		#	set ( OPTION_ON TRUE )
		#	set ( COMPILE_FLAG_OPTION "-DCMAKE_CXX_FLAGS:STRING=${COMPILE_FLAG}" )
		#endif ( "${COMPILE_FLAG}" STREQUAL "" )

		if    ( "${LINK_LIB}" STREQUAL "" )
			set ( EXTRA_OPTION "" )
		else  ( "${LINK_LIB}" STREQUAL "" )
			set ( OPTION_ON TRUE )
			set ( EXTRA_OPTION "-DLINK_LIBRARIES:STRING=${LINK_LIB}" )
		endif ( "${LINK_LIB}" STREQUAL "" )

		#MESSAGE ( STATUS "COMPILE_FLAG_OPTION ${COMPILE_FLAG_OPTION}" )
		#MESSAGE ( STATUS "EXTRA_OPTION ${EXTRA_OPTION}" )

		IF ( CROSS_COMPILING )

			if    ( OPTION_ON )
				try_compile ( ${RESULT_VAR} "${_bindir}" "${_SRCFILE}" CMAKE_FLAGS ${COMPILE_FLAG_OPTION} ${EXTRA_OPTION} )
			else  ( OPTION_ON )
				try_compile ( ${RESULT_VAR} "${_bindir}" "${_SRCFILE}" )
			endif ( OPTION_ON )

			IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )

				if    ( OPTION_ON )
					try_compile ( ${RESULT_VAR} "${_bindir}_fail" "${_SRCFILE_FAIL}" CMAKE_FLAGS ${COMPILE_FLAG_OPTION} ${EXTRA_OPTION} )
				else  ( OPTION_ON )
					try_compile ( ${RESULT_VAR} "${_bindir}_fail" "${_SRCFILE_FAIL}" )
				endif ( OPTION_ON )

			ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
		ELSE ( CROSS_COMPILING )
			if    ( OPTION_ON )
				try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}" "${_SRCFILE}" CMAKE_FLAGS ${COMPILE_FLAG_OPTION} ${EXTRA_OPTION} )
			else  ( OPTION_ON )
				try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}" "${_SRCFILE}" )
			endif ( OPTION_ON )

			IF ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
				SET ( ${RESULT_VAR} TRUE )
			ELSE ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
				SET ( ${RESULT_VAR} FALSE )
			ENDIF ( _COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR )
			IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )

				if    ( OPTION_ON )
					try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}_fail" "${_SRCFILE_FAIL}" CMAKE_FLAGS ${COMPILE_FLAG_OPTION} ${EXTRA_OPTION} )
				else  ( OPTION_ON )
					try_run ( _RUN_RESULT_VAR _COMPILE_RESULT_VAR "${_bindir}_fail" "${_SRCFILE_FAIL}" )
				endif ( OPTION_ON )

				IF ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
					SET ( ${RESULT_VAR} TRUE )
				ELSE ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
					SET ( ${RESULT_VAR} FALSE )
				ENDIF ( _COMPILE_RESULT_VAR AND _RUN_RESULT_VAR )
			ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL} )
		ENDIF ( CROSS_COMPILING )

		IF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE} )

			if    ( OPTION_ON )
				try_compile ( _TMP_RESULT "${_bindir}_fail_compile" "${_SRCFILE_FAIL_COMPILE}" CMAKE_FLAGS ${COMPILE_FLAG_OPTION} ${EXTRA_OPTION} )
			else  ( OPTION_ON )
				try_compile ( _TMP_RESULT "${_bindir}_fail_compile" "${_SRCFILE_FAIL_COMPILE}" )
			endif ( OPTION_ON )

			IF ( _TMP_RESULT )
				SET ( ${RESULT_VAR} FALSE )
			ELSE ( _TMP_RESULT )
				SET ( ${RESULT_VAR} TRUE )
			ENDIF ( _TMP_RESULT )
		ENDIF ( ${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE} )

		IF ( ${RESULT_VAR} )
			#MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME} -- works" )
			LIST ( APPEND CXX11_SUPPORTED_FEATURE_LIST ${FEATURE_NAME} )
		ELSE ( ${RESULT_VAR} )
			#MESSAGE ( STATUS "Checking C++11 support for ${_LOG_NAME} -- not supported" )
			LIST ( APPEND CXX11_UNSUPPORTED_FEATURE_LIST ${FEATURE_NAME} )
		ENDIF ( ${RESULT_VAR} )

		SET ( ${RESULT_VAR} ${${RESULT_VAR}} CACHE INTERNAL "C++11 support for ${_LOG_NAME}" )
	ENDIF ( NOT DEFINED ${RESULT_VAR} )

ENDMACRO ( checkFeature FEATURE_NAME FEATURE_NUMBER RESULT_VAR )

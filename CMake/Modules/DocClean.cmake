if ( TARGET docclean )
    # Target already available, do no create it then anymore
else ( TARGET docclean )
	add_custom_target ( docclean )
	add_custom_command (
		TARGET docclean
		COMMAND sh ${CMAKE_MODULE_PATH}/docclean.sh ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen
	)
endif ( TARGET docclean )
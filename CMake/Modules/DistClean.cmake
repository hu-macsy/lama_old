
add_custom_target ( distclean )
file ( GLOB_RECURSE BUILD_GLOB_RES ${CMAKE_BINARY_DIR}/* )
add_custom_command (
    TARGET distclean
    DEPENDS clean
    COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_GLOB_RES}
)
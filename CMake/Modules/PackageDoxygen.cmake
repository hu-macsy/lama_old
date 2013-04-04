### DOXYGEN DOCUMENTATION ###

if ( DOXYGEN_FOUND )

    set ( LAMA_DOC_DIR "${LAMA_SOURCE_DIR}/../doc" )
    set ( DOXYGEN_BUILD_ROOT "${LAMA_DOC_DIR}/doxygen" )
    filE ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command (
      OUTPUT ${DOXYGEN_BUILD_ROOT}/html/index.html
      COMMAND rm -rf ${DOXYGEN_BUILD_ROOT}
      COMMAND mkdir ${DOXYGEN_BUILD_ROOT}
      COMMAND ${DOXYGEN_EXECUTABLE} ${LAMA_DOC_DIR}/LAMA.Doxyfile
      DEPENDS ${LAMA_DOC_DIR}/LAMA.Doxyfile
      WORKING_DIRECTORY ${LAMA_DOC_DIR}
   )

   add_custom_target (
      doc
      DEPENDS
      ${DOXYGEN_BUILD_ROOT}/html/index.html
   )

   # Install the documentation generated at "make" time.

   # install ( DIRECTORY ${DOXYGEN_BUILD_ROOT}/ DESTINATION ${DOXYGEN_BUILD_ROOT}/html )

else ( DOXYGEN_FOUND )
    message ( WARNING "Not building system documentation because Doxygen not found." )
endif ( DOXYGEN_FOUND )
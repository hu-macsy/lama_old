###
 # @file PackageSphinx.cmake
 #
 # @license
 # Copyright (c) 2013
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
 # @brief configuration of sphinx
 # @author Jan Ecker
 # @date 15.05.2013
###

### DOXYGEN DOCUMENTATION ###

if ( SPHINX_FOUND )
    ### install ###
    set ( LAMA_DOC_DIR "${LAMA_SOURCE_DIR}/doc" )
    set ( SPHINX_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/doc/user" )
    file ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )
    
    configure_file( "${CMAKE_SOURCE_DIR}/doc/user/conf.py.in" "${CMAKE_SOURCE_DIR}/doc/user/conf.py" )
    configure_file( "${CMAKE_SOURCE_DIR}/doc/user/convert_json.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/user/convert_json.sh" )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command (
      OUTPUT ${SPHINX_BUILD_ROOT}/html/index.html
      COMMAND ${Sphinx-build_EXECUTABLE} -b html -d ${SPHINX_BUILD_ROOT}/doctrees user ${SPHINX_BUILD_ROOT}/html
      DEPENDS ${CMAKE_SOURCE_DIR}/doc/user/conf.py
      WORKING_DIRECTORY ${LAMA_DOC_DIR}
   )
   
   add_custom_command (
      OUTPUT ${SPHINX_BUILD_ROOT}/json/index.html
      COMMAND ${Sphinx-build_EXECUTABLE} -b json -d ${SPHINX_BUILD_ROOT}/doctrees user ${SPHINX_BUILD_ROOT}/json
      COMMAND chmod +x ${SPHINX_BUILD_ROOT}/convert_json.sh
      COMMAND ${SPHINX_BUILD_ROOT}/convert_json.sh
      DEPENDS ${CMAKE_SOURCE_DIR}/doc/user/conf.py
      DEPENDS ${SPHINX_BUILD_ROOT}/convert_json.sh
      WORKING_DIRECTORY ${LAMA_DOC_DIR}
   )
   
   #$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
   # cd $(BUILDDIR)/json; \
   # echo `pwd`; \
   # find ./ -type f -exec sed -i 's/\.\.\/\_images/fileadmin\/LAMA\/json\/_images/g' {} \; \
   # && find ./ -type f -exec sed -i 's/<img/<image/g' {} \; \
   # && find ./ -type f -exec sed -i 's/\.\.\/fileadmin/fileadmin/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"configuration\/\\\"/\\\"configuration\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"installation\/\\\"/\\\"installation\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"tutorial\/\\\"/\\\"tutorial\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"lecture\/\\\"/\\\"lecture\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"reference\/\\\"/\\\"reference\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"testing\/\\\"/\\\"testing\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"benchmarks\/\\\"/\\\"benchmarks\/index\/\\\"/g' {} \; \
#    && find ./ -type f -exec sed -i 's/\\\"developer\/\\\"/\\\"developer\/index\/\\\"/g' {} \;\
#    && find ./ -type f -exec sed -i 's/\\\"solver\/\\\"/\\\"solver\/index\/\\\"/g' {} \;
#    @echo
#    @echo "Build finished; now you can process the JSON files."

   add_custom_target (
      userdoc
      DEPENDS
      ${SPHINX_BUILD_ROOT}/html/index.html
   )
   add_custom_target (
      userdoc_json
      DEPENDS
      ${SPHINX_BUILD_ROOT}/json/index.html
   )

   # Install the documentation generated at "make" time.

   # install ( DIRECTORY ${DOXYGEN_BUILD_ROOT}/ DESTINATION ${DOXYGEN_BUILD_ROOT}/html )

else ( SPHINX_FOUND )
    if ( LAMA_CMAKE_VERBOSE )
        message ( STATUS "Not building user documentation because Sphinx not found." )
    endif ( LAMA_CMAKE_VERBOSE )
endif ( SPHINX_FOUND )
### set include directories ###

include_directories( . yaml/ )

### add subdirectories ###

add_subdirectory( ./src )

file ( COPY BenchmarkRunner DESTINATION ${CMAKE_BINARY_DIR}/framework )
file ( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/framework/BenchmarkFramework/.log )


### set source files ###

set ( PYTHON_FILES
       __init__.py
       Benchmark.py
       BenchmarkOut.py
       Column.py
       Configure.py
       Enum.py
       Exceptions.py
       InputSet.py
       LaunchConfiguration.py
       Options.py
       Run.py
       Sort.py
       texttable.py
     )
     
#
# export file: copy it to the build tree on every build invocation and add rule for installation
# solution from:
# http://stackoverflow.com/questions/697560/cmake-how-to-copy-directory-from-source-tree-to-binary-tree
#
function ( cm_export_file FILE DEST )
  if ( NOT TARGET export-files )
    add_custom_target ( export-files ALL COMMENT "Exporting files into build tree" )
  endif ( NOT TARGET export-files )
  get_filename_component ( FILENAME "${FILE}" NAME )
  add_custom_command ( TARGET export-files COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${DEST}/${FILENAME}")
  install ( FILES "${FILE}" DESTINATION "${DEST}" )
endfunction ( cm_export_file )

foreach( PYTHON_FILE ${PYTHON_FILES} )
    cm_export_file ( BenchmarkFramework/${PYTHON_FILE} BenchmarkFramework/ )
endforeach ( PYTHON_FILE ${PYTHON_FILES} )

install ( FILES "BenchmarkRunner" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin )

file ( COPY BenchmarkFramework/yaml DESTINATION ${CMAKE_BINARY_DIR}/framework )

# Apply export/install command also to all yaml files

set ( YAML_DIR ${CMAKE_SOURCE_DIR}/framework/BenchmarkFramework/yaml )
file ( GLOB YAML_FILES RELATIVE ${YAML_DIR} ${YAML_DIR}/*.py )
foreach( PYTHON_FILE ${YAML_FILES} )
    cm_export_file ( BenchmarkFramework/yaml/${PYTHON_FILE} BenchmarkFramework/yaml/ )
endforeach ( PYTHON_FILE ${YAML_FILES} )

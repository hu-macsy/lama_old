
set (TESTSUPPORT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include PARENT_SCOPE)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_summary.xslt  DESTINATION ${CMAKE_BINARY_DIR}/testsupport )
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py  DESTINATION ${CMAKE_BINARY_DIR}/testsupport )
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/testsupport.py  DESTINATION ${CMAKE_BINARY_DIR}/testsupport )

add_library(testsupport STATIC src/GlobalTempDir.cpp src/uniquePath.cpp)
target_include_directories(testsupport PRIVATE include)

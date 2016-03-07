BOOST
^^^^^

`Boost`_ provides free portable C++ source libraries.

.. _Boost: http://www.boost.org

LAMA uses some Boost libraries and therefore you must have an actual version (1.34 or later) installed on your machine.
The amount of used Boost libraries dependens on the compiler. If the compiler provides support for C++11 fewer parts of Boost
will be used. 

The following only header libraries are used by LAMA (obsolete by using a C++11 compiler):

  - **Smart pointers**: avoid to keep track of ownership of dynamically allocated memory
  - **Function**: object wrappers for deferred calls or callbacks.
  - **Bind**: allows to bind arguments for new function pointers
  - **Boost Spirit Parser Framework** (object oriented recursive descent parser generator) (optional) 

Beside some only header libraries, the following *compiled* libraries are used by LAMA:

  - **thread** (portable C++ multi-threading, mandatory)
  - **unit_test_framework** and **regex** (program and full unit testing, recommended)

.. - **program_options** (easy access to options of a program call, optional)

Many linux installations provide an actual release of Boost and if Boost is installed, the LAMA configuration should
have no problems to find it.

If for any reasons no actual Boost installation is available, you must download and install it. 
Please make sure that you build also the dynamic library versions. After installation of Boost you can tell cmake 
the location of installation by the variable BOOST_ROOT or by setting an environment variable::

    export BOOST_ROOT=<path/to/boost/installation>
    cmake ..

or::

    cmake -DBOOST_ROOT=<path/to/boost/installation> ...

Via ccmake you can verify that the Boost variables needed for LAMA have correct values::

    BOOST_ROOT                        <path/to/boost/installation>
    Boost_INCLUDE_DIR                 <path/to/boost/installation>/include
    Boost_LIBRARY_DIRS                <path/to/boost/installation>/lib
    Boost_REGEX_LIBRARY               <path/to/boost/installation>/lib/libboost_regex.so
    Boost_SYSTEM_LIBRARY              <path/to/boost/installation>/lib/libboost_system.so
    Boost_THREAD_LIBRARY              <path/to/boost/installation>/lib/libboost_thread.so
    Boost_UNIT_TEST_FRAMEWORK_LIBRARY <path/to/boost/installation>/lib/libboost_unit_test_framework.so
    
..    Boost_PROGRAM_OPTIONS_LIBRARY     <path/to/boost/installation>/lib/libboost_program_options.so


There are known issues with some older Boost Installations with there own Boost.cmake definition.
If you have error message looking like::

    CMake Error at /usr/lib64/boost/Boost.cmake:536 (message):
    The imported target "boost_date_time-static" references the file

      "/usr/lib64/lib64/libboost_date_time.a"

    but this file does not exist.  Possible reasons include:

    * The file was deleted, renamed, or moved to another location.

    * An install or uninstall procedure did not complete successfully.

    * The installation package was faulty and contained

     "/usr/lib64/boost/Boost.cmake"

    but not all the files it references.

try using -DBoost_NO_BOOST_CMAKE=TRUE
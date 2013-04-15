.. _section_boost:

BOOST
-----

`Boost`_ provides free portable C++ source libraries.

.. _Boost: http://www.boost.org

LAMA uses some Boost libraries and therefore you must have an actual version (1.34 or later) installed on your machine.

The following only header libraries are used by LAMA:

  - **Smart pointers**: avoid to keep track of ownership of dynamically allocated memory
  - **Function**: object wrappers for deferred calls or callbacks.
  - **Bind**: allows to bind arguments for new function pointers
  - **Boost Spirit Parser Framework** (object oriented recursive descent parser generator) (optional) 

Beside some only header libraries, the following *compiled* libraries are used by LAMA:

  - **thread** (portable C++ multi-threading, mandatory)
  - **unit_test_framework** and **regex** (program and full unit testing, recommended)
  - **program_options** (easy access to options of a program call, optional)

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

    BOOST_ROOT                     <path/to/boost/installation>
    Boost_INCLUDE_DIR              <path/to/boost/installation>/include
    Boost_LIBRARY_DIRS             <path/to/boost/installation>/lib
    Boost_PROGRAM_OPTIONS_LIBRARY  <path/to/boost/installation>/lib/libboost_program_options.so
    Boost_THREAD_LIBRARY           <path/to/boost/installation>/lib/libboost_thread.so
    Boost_UNIT_TEST_FRAMEWORK_LIBR <path/to/boost/installation>/lib/libboost_unit_test_framework.so


.. _section_boost:

BOOST
-----

`Boost`_ provides free portable C++ source libraries.

.. _Boost: http://www.boost.org

LAMA uses some Boost libraries and therefore you must have an actual version (1.34 or later) installed on your machine.

The following only header libraries are used by LAMA:

  - Smart pointers: avoid to keep track of ownership of dynamically allocated memory
  - Function: object wrappers for deferred calls or callbacks.
  - Bind: allows to bind arguments for new function pointers
  - Boost Spirit Parser Framework (object oriented recursive descent parser generator) 

Beside some only header libraries, the following 'compiled' libraries are used by LAMA:

  - Thread (portable C++ multi-threading, mandatory)
  - Test (program and full unit testing, recommended)
  - Program Options (easy access to options of a program call, optional)

Many linux installations provide an actual release of Boost and if Boost is installed, 
the LAMA configuration should have no problems to find it.

If for any reasons no actual Boost installation is available, you must download and install it. 
Please make sure that you build also the dynamic library versions. After installation of Boost you can tell cmake 
the location of installation by the variable BOOST_ROOT as described above or by setting an environment variable::

    export BOOST_ROOT=<path-to-boost-installation>
    cmake ..

or::

    cmake -DBOOST_ROOT=<path-to-boost-installation> ...

Via ccmake you can verify that the Boost variables needed for LAMA have correct values::

    BOOST_ROOT                     /home/brandes/local/boost_1_46_0
    Boost_INCLUDE_DIR              /home/brandes/local/boost_1_46_0/include
    Boost_LIBRARY_DIRS             /home/brandes/local/boost_1_46_0/lib
    Boost_PROGRAM_OPTIONS_LIBRARY  /home/brandes/local/boost_1_46_0/lib/libboost_program_options.so
    Boost_THREAD_LIBRARY           /home/brandes/local/boost_1_46_0/lib/libboost_thread.so
    Boost_UNIT_TEST_FRAMEWORK_LIBR /home/brandes/local/boost_1_46_0/lib/libboost_unit_test_framework.so


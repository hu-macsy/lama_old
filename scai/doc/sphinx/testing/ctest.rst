Using CTest
-----------

CTest is part of Cmake and helps to integrate various testing executable with various runtime configurations to the existing buildsystem.

In an existing CMakeLists.txt file there are two Makros to use:

To activate a testing facility of CMake/CTest you have to use :

.. code-block:: bash

	ENABLE_TESTING()

After this Makro you are able to add different test executables by using:

.. code-block:: bash

	ADD_TEST(<NAME> <COMMAND>)

e.g.:

.. code-block:: bash

	ADD_TEST( TestRunName ./testrun )

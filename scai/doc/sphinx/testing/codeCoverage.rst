Code Coverage on LAMA
---------------------

This article shows how to measure code coverage on LAMA. Using Lcov is a good option to do this measurement
and to create a html-structure as an output.

How to measure code coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requirements
^^^^^^^^^^^^

To be able to create gcno and gcda files, it is required to compile a project with the flags "--coverage" or
"-fprofile-arcs -ftest-coverage". 

How to
^^^^^^ 

Executing the following commands in your build directory and the code coverage will be measured:

.. code-block:: bash

	lcov --directory .. --zerocounters

Now run your program.

.. code-block:: bash
	
	lcov --directory .. --capture --output-file=data.info
	genhtml data.info


Usability in LAMA
^^^^^^^^^^^^^^^^^

LAMA itself offers the environment variable CODE_COVERAGE=ON to enable all needed flags.
While building LAMA with this option, a script cca.sh will be copied into the testing directory. This script invokes 
all neccessary steps to measure the coverage. The result will be a HTML-structure showing the analysis of the 
measurement.

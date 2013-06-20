Continuous Integration in LAMA
==============================

Continuous Integration (CI) provides different advantages for a software development project: Centralisation of all
necessary data of this project, automatization of processes, like nightly builds or summing up testresults, an early
elimination of errors and propagation of different information.

We have decided to use Jenkins as a CI-Server in LAMA. To obtain all needed facilities, we extend this server with some
plugins.

- Hudson Post Build task : Let you execute some shell commands after the buildprocess

- xUnit Plugin : This plugin gives the possibility to collect and sum up testresults of xml-documents created by the
  used test framework.
  
- Hudson Cmake Builder : Let you build automatically a cmake-based project

Configuration of a LAMA-Job
---------------------------

General option:

Source-Code-Management : SVN  with repository url

::

	svn+ssh://**USERNAME**@svn.code.sf.net/p/libama/code/trunk

It is recommended to activate the option "Poll SCM" as a build-trigger. The schedule should be configured aswell to
check the SCM just e.g. every 5 minutes, " \*/5 * * * * ". This reduces the number of created builds, if many commits
appear at the same time.

The Cmake Build plugin contains the following values:

- Source directory : src/
- Build directory : build/

Advanced:

- Other CMake Arguments (as an example) : -DBOOST_ROOT=/usr/include/boost -DLAMA_USE_CUDA=OFF -DLAMA_USE_OPENCL=OFF

After building a project, all tests should be executed to varify the correctness of the calculations.
This can be done by writing a shell script in a post build task:

::

	cd build/test
	echo "Exception = ERROR" >> config
	export LAMA_LOG=config
	export LAMA_UNSUPPORTED=IGNORE

	./lama_test --log_level=test_suite 

Prospective option
------------------

Boost.Test version 1.42 and higher offers the runtime parameter log\_sink=file to push the output into a specific file. 
This option can be used to push the XML-output created by the runtime parameter log\_format=XML into a \*.xml file, that 
can be read from the xUnit-Plugin. This plugin creates a diagram on the mainpage of the LAMA-Job to show the working and
failing tests.

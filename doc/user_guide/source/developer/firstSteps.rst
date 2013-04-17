First Steps
===========

This description should give you a good start in programming with and in LAMA. It will show you the requirements for a 
convenient implementation process, how to work with the versioning tool and how to properly build and test LAMA.

Requirements
------------

  - programming environment: Eclipse IDE C/C++
  - versioning tool: Git
  - build tool: CMake
  - test framework: Boost test
  - documentation tool: Doxygen (system docu) and Sphinx (user docu --> web page)

Eclipse
^^^^^^^
We are using the Eclipse IDE for C/C++ as our standard programming environment. It is well known, so most programmers
are comfortable with it and it offers a couple of nice features as there are plugins for versioning, cmake and ReST
(Sphinx) editor as well as a build-in debugger or code-formatting support.  


- Java SDK and Eclipse IDE C/C++ have to be both either 32bit or 64bit (esp. Linux)! 
- Please, also refer to the Eclipse documentation in order to avoid problems with non-suitable java versions
- The environment variable **JAVA_HOME** has to exist and point to the right location

.. code-block:: bash

	JAVA_HOME=<install/dir/java>/bin
	
We recommend Eclipse CDT `Indigo`_ now at least with the `EGit`_ plugin. Optional you can install the
`CMake <http://www.cthing.com>`_ and `ReST`_ editor. 

.. _Indigo: http://www.eclipse.org/downloads/packages/eclipse-ide-cc-developers-includes-incubating-components/indigosr2
.. _EGit: http://www.eclipse.org/egit/
.. _ReST: http://resteditor.sourceforge.net/

GIT
^^^

Our versioning tool of choice is `Git`_. The repository is hosted at sourceforge.
You can check out the source code via Git where the HTTP access is as follows:

.. _Git: http://git-scm.com/

.. code-block:: bash

    git clone http://git.code.sf.net/p/libama/git lama

Code developers must register at `Sourceforge`_ and join the LAMA developer team.
Place your rsa-key at Account --> Services to be able for checking in.
You should download the software as follows:

.. _Sourceforge: http://sourceforge.net/

.. code-block:: bash

   git clone ssh://<your_user_id>@git.code.sf.net/p/libama/git lama

Using the EGit Eclipse plugin just import a Git project and use ``ssh://<your_user_id>@git.code.sf.net/p/libama/git`` as
repository location.

CMake
^^^^^

We are using `CMake <http://www.cmake.org/>`_ for creating our Makefiles (Linux) and Visual Studio project (Windows).
Therefore CMake (in Version 2.8 or greater) have to be installed and be present in the standard executable path.

Doxygen
^^^^^^^

For our system documentation we make use of `Doxygen`_. So if API documentation is supposed to be generated Doxygen has
to be present on the system. To build it just call ``make doc`` in the build directory after the configuration step.

.. _Doxygen: http://www.doxygen.org

Sphinx
^^^^^^

`Sphinx`_ is used for our user documentation. Sphinx offers a wide support of different targets. So you can create a
local webbased documentation by calling ``make html`` or a pdf by ``make pdflatex`` in **LAMA_ROOT/doc/user_guide/**.
The Sphinx source are located in LAMA_ROOT/doc/user_guide/source/, refered C++-sources of the lecture and the tutorial
you can find in the subfolders of LAMA_ROOT/doc/user_guide/cpp_source/.

.. _Sphinx: http://sphinx-doc.org/
.. _main-page_logging:

############
SCAI Logging
############

***********
Description
***********

Inserting log statements into code is a low-tech method for debugging it. It may also be the only way because
debuggers are not always available or applicable. This is usually the case for multithreaded applications and
distributed applications at large.

Experience indicates that logging is an important component of the development cycle. It offers several
advantages. It provides precise context about a run of the application. Once inserted into the code, the
generation of logging output requires no human intervention. Moreover, log output can be saved in persistent
medium to be studied at a later time. In addition to its use in the development cycle, a sufficiently rich
logging package can also be viewed as an auditing tool.

Logging does have its drawbacks. It can slow down an application. If too verbose, it can cause scrolling
blindness. To alleviate these concerns, our logging module is designed to be reliable, fast and extensible.
Since logging is rarely the main focus of an application, the API strives to be simple to understand and to
use.

* Contains macros for defining loggers and logging calls

We decided to use an own simple logging facilitate to decrease dependencies with other software packages.
The SCAI logging library is a very convenient and efficient library that supports logging in C++ applications.
It has very low overhead and has quite powerful mechanisms for generating logging output.

It is used in all other SCAI projects but can also be used in other applications.

Currently, logging is always done on standard output. It is not possible to write logging messages in different 
output streams.

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   Integration
   Compilation
   Runtime

*******
Example
*******

.. toctree::
   :titlesonly:
   
   Example

************
Dependencies
************

Internal dependencies:

* :ref:`SCAI Common<scaicommon:main-page_common>` provides the 
  :ref:`Thread class<scaicommon:Thread>` for naming of threads

*****
Usage
*****

A CMake module is provided that finds the include dir and the logging library and sets the compile flag
for the desired logging level.

* compile flag `-DSCAI_LOG_LEVEL_xxx` must be set at compile time to set the logging level
* environment variable ``SCAI_LOG`` specifies either the logging level or a configuration file
  
************
Related Work
************

The following three versions of logging libraries are used in most other C++ projects:

- |log4cpp|
- |log4cxx|
- |boostlog|

.. |log4cpp| raw:: html

   <a href="http://log4cpp.sourceforge.net" target="_blank"> log4cpp </a>

.. |log4cxx| raw:: html

   <a href="http://logging.apache.org/log4cxx" target="_blank"> log4cxx </a>

.. |boostlog| raw:: html

   <a href="http://boost-log.sourceforge.net/libs/log/doc/html/index.html" target="_blank"> Boost logging </a>

There are certain advantages and disadvantges for each of the libraries.
There is a comparative survey `here`__ for the different facilities but this survey does not give any hints
about performance, usability, and availability. Here are some impressions for each of the libraries:

__ http://log4cpp.hora-obscura.de/index.php/LoggingLibraryForCpp

- log4cpp is very simple, easy to install, but not very performant.

- log4cxx is the most popular, very performant, highest functionality, but requires the Apache runtime system
  and utilities to be installed.
  
- Boost logging might be integrated in one of the next versions of Boost but not yet; it has only include
  files but no library itself; runtime configuration files are not supported.

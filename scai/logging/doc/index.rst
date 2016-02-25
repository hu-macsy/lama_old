.. _main-page_logging:

Introduction
============

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

Specification
-------------

* Contains macros for defining loggers and logging calls
* Internal dependencies: common

SCAI Logging
------------

We decided to use an own simple logging facilitate to decrease dependencies with other software packages.
The macros make it possible to implement the logging with one arbitrary library and may be later
with an own version.




.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   Integration
   Compilation
   Runtime
   Example
   RelatedWork


Summary
-------

The SCAI logging library is a very convenient and efficient library that supports logging in C++ applications.
It has very low overhead and has quite powerful mechanisms for generating logging output.

It is used in all other SCAI projects coming with the distribution but can also be used in other applications.
A CMake module is provided that finds the include dir and the logging library and sets the compile flag
for the desired logging level.

Currently, logging is always done on standard output. It is not possible to write logging messages in different 
output streams.



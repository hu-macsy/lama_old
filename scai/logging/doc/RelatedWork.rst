Related Work
============

The following three versions of logging libraries are used in most other C++ projects:

- `log4cpp`_
- `log4cxx`_
- `boost logging <http://boost-log.sourceforge.net/libs/log/doc/html/index.html>`_

.. _log4cpp: http://log4cpp.sourceforge.net/
.. _log4cxx: http://logging.apache.org/log4cxx/

There are certain advantages and disadvantges for each of the libraries.
There is a comparative survey `here`__ for the different facilities but this survey does not give any hints
about performance, usability, and availability. Here are some impressions for each of the libraries:

__ http://log4cpp.hora-obscura.de/index.php/LoggingLibraryForCpp

- log4cpp is very simple, easy to install, but not very performant.

- log4cxx is the most popular, very performant, highest functionality, but requires the Apache runtime system
  and utilities to be installed.
  
- Boost logging might be integrated in one of the next versions of Boost but not yet; it has only include
  files but no library itself; runtime configuration files are not supported.
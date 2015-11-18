.. _main-page:

SCAI Common
===========

* Provides different functionality needed in nearly all SCAI projects
* Provides coding conventions for using smart pointers
* Provides coding conventions for exception handling
* Provides macros for assertions
* Provides macros to deal with different platforms
* Own classes for platform-specific stuff (e.g. threads, ... )

Classes
-------

====================    ==========================================
Class                   Description
====================    ==========================================
Exception               Error handling, call stack
Factory                 Template class for Factory
Factory1                Factory, but create with additional argument
Thread                  Basic stuff to deal with multithreading (uses pThreads)
ScalarType              Enumeration type for supported value types (allows registration for different types in factories, interfaces, ...)
Walltime                Simple and efficient walltime measuring
Printable               Base class to support stream output
NonCopyable             Disable default copy constructors
shared_ptr              Either boost::shared_ptr or std::shared_ptr
function, bind          Eiter boost::function or std::function
LibModule               Load/Unload of Library Modules (dynamic libraries)
ContextType             Enum types ContextType and AccessKind
Settings                Access to environment variables
OpenMP                  Dummy routines if OpenMP is disabled
SCAIType
====================    ==========================================
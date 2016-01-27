.. _main-page:

SCAI Common
===========

The common library contains different functionality needed in nearly all SCAI projects.
These are some utitily classes (e.g. for exceptions, factories, timing, ...), some
classes to abstract from platform-specific features (e.g. threads, loading library modules, ...),
and some kind of coding conventions (e.g. for using smart pointers, function types, ...).
Beside these more general concepts that might also be used for other libraries, the common
library provides some LAMA specific enumeration types (supported value types, different context
types, ...).

All classes and type defintions are done in the namespace ``scai/common``.

SCAI Types
----------

====================    ==========================================
Class                   Description
====================    ==========================================
ContextType             Enumeration types ContextType and AccessKind
ScalarType              Enumeration type for supported value types (allows registration for different types in factories, interfaces, ...)
TypeTrait               Struct with all specific stuff for any supported arithmetic value type in matrix/vector operations
SCAIType                Supported arithmetic types for template instantiations
====================    ==========================================

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   Types
   TypeTrait
   Constants
   SCAITypes
   Complex

Common Classes
--------------

====================    ==========================================
Class                   Description
====================    ==========================================
Exception               Error handling, call stack
Factory                 Template class for Factory
Factory1                Factory, but create with additional argument
Thread                  Basic stuff to deal with multithreading (uses pThreads)
Walltime                Simple and efficient walltime measuring
Printable               Base class to support stream output
NonCopyable             Disable default copy constructors
shared_ptr              Either boost::shared_ptr or std::shared_ptr
function, bind          Either boost::function or std::function
LibModule               Load/Unload of Library Modules (dynamic libraries)
ContextType             Enum types ContextType and AccessKind
Settings                Access to environment variables
OpenMP                  Dummy routines if OpenMP is disabled
====================    ==========================================

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   Exception
   Printable
   Factory
   Timing
   LibModule
   Thread
   Settings
   SmartPointers
   Function
   OpenMP

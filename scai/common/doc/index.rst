.. _main-page_common:

###########
SCAI Common
###########

The common library contains different functionality needed in nearly all SCAI projects.
These are some utitily classes (e.g. for exceptions, factories, timing, ...), some
classes to abstract from platform-specific features (e.g. threads, loading library modules, ...),
and some kind of coding conventions (e.g. for using smart pointers, function types, ...).
Beside these more general concepts that might also be used for other libraries, the common
library provides some LAMA specific enumeration types (supported value types, different context
types, ...).

All classes and type defintions are done in the namespace ``scai/common``.

*****************
Enumeration Types
*****************

Within all SCAI projects, enum types are used when the number of possible values is limited and when
the values have a special meaning.

The following enum types are defined in common as they are use within multiple projects.

====================    ==========================================
Class                   Description
====================    ==========================================
:ref:`ContextType`      Enumeration types ContextType and AccessKind
:ref:`ScalarType`       Enumeration type for supported value types (allows registration for different types in factories, interfaces, ...)
:ref:`ReductionOp`      Enumeration type for different binary operators used in reductions
====================    ==========================================

****************
Arithmetic Types
****************

LAMA is a library for mathematical operations and here in common it is defined for which arithmetic types
templated classes and operations will be instantiated, where the template arguments stand for an arithmetic type.
An arithmetic type that is used for the instantation of these classes must provide a certain number of 
operations, e.g. +, -, \*, /, and so on. Furthermore some mathematical operations and type properties must be 
provided. Some additional structures are used to provide these functionalities in a consistent way. 

====================    ==========================================
Class                   Description
====================    ==========================================
:ref:`SCAITypes`        Supported arithmetic types for template instantiations
:ref:`TypeTrait`        Struct with all specific stuff for any supported arithmetic value type in matrix/vector operations
:ref:`Math`             Wrapper for mathematical operations (like those in cmath) to use them in templated code
:ref:`Complex`          Complex numbers which cannot only be used on host, but also on CUDA and MIC devices.
:ref:`Constants`        Operations to compare value to a machine specific eps
====================    ==========================================

**************
Common Classes
**************

====================         ==========================================
Class                        Description
====================         ==========================================
:ref:`Exception`             Error handling, call stack
:ref:`Factory`               Template class for Factory
:ref:`Factory1`              Factory, but create with additional argument
:ref:`Thread`                Basic stuff to deal with multithreading (uses pThreads)
:ref:`Walltime`              Simple and efficient walltime measuring
:ref:`Printable`             Base class to support stream output
:ref:`NonCopyable`           Disable default copy constructors
:ref:`LibModule`             Load/Unload of Library Modules (dynamic libraries)
====================         ==========================================

***************
Common Concepts
***************

The following stuff stands for general concepts that is or might be used in all SCAI projects.

====================         ==========================================
Name                         Description
====================         ==========================================
:ref:`SmartPointers`         Smart pointers that deallocate object with destructor: unique_ptr or shared_ptr or unique_ptr)
:ref:`Assertion`             Assertion checking, which can be compiled out of code
:ref:`Function`              Function objects that might also be created with bound arguments
:ref:`Settings`              Access to environment variables
:ref:`OpenMP`                Dummy routines if OpenMP is disabled
====================         ==========================================

***********************
Common Classes for CUDA
***********************

Some general stuff used for CUDA is also part of the common project.

.. toctree::
   :titlesonly:
   :maxdepth: 2

   CUDAError
   CUDACtx
   CUDAAccess

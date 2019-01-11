.. _main-page_common:

###########
SCAI Common
###########

***********
Description
***********

The common library contains different functionality needed in nearly all SCAI modules.
These are some utitily classes (e.g. for exceptions, factories, timing, ...), some
classes to abstract from platform-specific features (loading library modules, ...)
Beside these more general concepts that might also be used for other libraries, the common
library provides some LAMA specific enumeration types (supported value types, different context
types, ...).

Nearly all classes and type defintions are set in the namespace ``scai/common``, only
the most important type definitions used for indexes and values in matrices are set
in the namespace ``scai``.

************************
Common Library Reference
************************

Utility Classes
---------------

The following classes provide some utitilites used in nearly all other
SCAI projects.

====================         ==========================================
Class                        Description
====================         ==========================================
:ref:`Factory`               Template class for Factory
:ref:`Thread`                Giving threads a unique name for identification
:ref:`Walltime`              Simple and efficient walltime measuring
:ref:`Grid`                  Definition of a rectangular grid
:ref:`Stencil`               Definition of a stencil
:ref:`Printable`             Base class to support stream output
:ref:`NonCopyable`           Disable default copy constructors
:ref:`LibModule`             Load/Unload of Library Modules (dynamic libraries)
====================         ==========================================

.. toctree::
   :hidden:

   Factory
   Thread
   Walltime
   Grid
   Stencil
   Printable
   NonCopyable
   LibModule

Common Concepts
---------------

The following stuff stands for general concepts that are or might be used in all SCAI projects.

====================         ==========================================
Name                         Description
====================         ==========================================
:ref:`Exception`             Error handling, call stack
:ref:`Assertion`             Assertion checking, which can be compiled out of code
:ref:`Settings`              Access to environment variables
:ref:`OpenMP`                Dummy routines if OpenMP is disabled
:ref:`TypeList`              Meta programming schemes to deal with list of types
====================         ==========================================

.. toctree::
   :hidden:

   Exception
   Assertion
   SmartPointers
   Settings
   OpenMP
   TypeList

Macros
------

The common project provides a lot of useful macros. Usually, each macro is
defined within a corresponing header file in the subdirectory ``macros``.

.. toctree::
   :titlesonly:
   :maxdepth: 2

   Macros

Enumeration Types
-----------------

Within all SCAI projects, enum types are used when the number of possible values is limited and when
the values have a special meaning.

The following enum types are defined in common as they are use within multiple projects.

====================    ==========================================
Class                   Description
====================    ==========================================
:ref:`ContextType`      Enumeration type ContextType for supported devices
:ref:`AccessKind`       Enumeration type AccessKind for read/write accesses
:ref:`ScalarType`       Enumeration type for supported value types (allows registration for different types in factories, interfaces, ...)
====================    ==========================================

.. toctree::
   :hidden:

   ContextType
   AccessKind
   ScalarType

Arithmetic Types
----------------

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
:ref:`Complex`          Complex numbers which cannot only be used on host, but also on CUDA.
:ref:`Constants`        Operations to compare value to a machine specific eps
:ref:`BinaryOp`         Enumeration type for binary operations
:ref:`UnaryOp`          Enumeration type for unary operations
:ref:`CompareOp`        Enumeration type for comparison operations
====================    ==========================================

.. toctree::
   :hidden:

   SCAITypes
   TypeTrait
   Math
   Complex
   Constants
   BinaryOp
   UnaryOp
   CompareOp

Common Classes for CUDA
-----------------------

Some general stuff used for CUDA is also part of the common project.

====================         ==========================================
Name                         Description
====================         ==========================================
:ref:`CUDACtx`               CUDA context
:ref:`CUDAAccess`            Access to CUDA context
:ref:`CUDAError`             Error handling for CUDA
====================         ==========================================

.. toctree::
   :hidden:

   CUDAError
   CUDACtx
   CUDAAccess

*****
Usage
*****

* Compile flag must be set ``SCAI_ASSERT_LEVEL_DEBUG``
* When the native C++ compiler does not support the C+11 standard,
  the Boost header files are needed and the include path
  for the Boost header files must be specified for the compilation.

Using the Settings class requires that the command line arguments are parsed.

.. code-block:: c++

    int main( int argc, const char* argv[] )
    {
        common::Settings::parseArgs( argc, argv );
        ...
    }

The following environment variables are used by the COMMON library:

* ``SCAI_UNSUPPORTED`` (ignore, warn, or error) when using the macro ``SCAI_UNSUPPORTED``

************
Dependencies
************

The common project is on the lowest level of the SCAI project hierarchy.
Therefore it does not depend on any other SCAI project.

These are the external dependencies:

* When the C++11 standard is not supported, header libraries of :ref:`Boost` must be available
* When the C++11 standard is not supported, a :ref:`PThread` library must be available
  for the implementation of the Thread class.
* :ref:`CUDA` toolkit 

.. toctree::
   :hidden:

   Boost
   PThread
   CUDA

************
Related Work
************

* `Boost Libraries <http://www.boost.org>`_, some functionality has been taken over in the C++11 standard.
* Macros and Meta Programming are techniques that are well described in
  "C++ Template Metaprogramming", by David Abrahams and Aleksey Gurtovoy. Copyright (c) 2005 by Pearson Education, Inc. 
* The TypeList concept is well documented in "Modern C++ Design: Generic Programming and Design Patterns Applied",
  by Andrei Alexandrescu. Copyright (c) 2001 by Addison Wesley
* Factory stuff as dynamic extension, similiar to module conecpt of Python
* The Curiously recurring template pattern as an idiom of C++ is well described
  `here <https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_


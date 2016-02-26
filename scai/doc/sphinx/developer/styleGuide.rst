.. _styleGuide:

Style Guide
-----------

This style guide should improve the readability and help to maintain a common style of the LAMA source code. In general
it should be considered as a guide line not strict rule set. To fulfill these requirements and to make it easier to
flow each rule should have a motivation. This style guide is mainly intended for C++ therefore some exceptions may
apply to source code written in other languages.

General Remarks
^^^^^^^^^^^^^^^

- All source files should compile with all enabled warnings without emitting any warning message.
  This rule also applies for the generation of the documentation with doxygen and CMake. Rationale:
  compiler warnings in general are meaningful and can be hints for bugs that might occur at runtime,
  so they should be taken into account. Even if a single compiler warning can be savely ignored it
  should be eliminated anyhow because otherwise important warning can be overlooked within a growing
  number of unimportant warnings.

- The Unit test should run warning free with valgrind. To suppress warnings from auxiliary libraries
  a suppressions file is provided (see :ref:`tools`). Please keep in mind that memory
  leaks and access violations are simply bugs that need to be fixed. Occasionally occurring warnings
  that are no bugs should be fixed also, because important warnings could be overlooked and they make
  it harder for LAMA users to utilize a tool like valgrind for their own purposes.

- The project should be complied regularly (at least with each release) with different compilers, of
  course warning free. This rule keeps the source code portable and helps to detect possible bugs
  earlier because different compilers emit different warnings. This rule implies that no compiler
  specific language extensions should be used and that all source code should be standard compliant.

- Whereever it is reasonable and possible we insert assertions. The rationale for this rule is
  that bugs are easier to fix if they are detected earlier, both in terms of program execution and
  development time.

Formatting
^^^^^^^^^^

The rationale for the formatting rules is code readability and the maintenance of a common style within LAMA, if no
other rationale is given.

- Only one statement is allowed per source line, e.g. **++i;++j;** is not allowed. This rule reduces the posibility to
  overlook a statement.

- The length of a row should not exceed 120 characters. Methods with long argument list should be
  formatted in the following way

::
	
   void veryLongMethodName(
       ValueType longArgumentName1,
       ValueType longArgumentName2 );
    	
or this way

::

   template<typename T, typename Allocator>
       inline const Expression<
            Expression< T,
                CUDAVector<T, Allocator>,
                    Times>,
            Expression< T,
                Expression<
                    CUDACSRSparseMatrix<T, Allocator>,
                        CUDAVector<T, Allocator>,
                            Times>,
                        Times>,
                    Plus>
        operator+(
             const CUDAVector<T, Allocator>& arg1,
             const Expression<
                 CUDACSRSparseMatrix<T, Allocator>,
                 CUDAVector<T, Allocator>,
                 Times>& arg2 );
                  
This rule should inhibit the necessarity to scroll horizontally within source code editors and
allow the usage of different editors also on the command line (e.g. **vi**)

- Indention is done with four whitespace instead of a tabulator. This makes the source code more
  readable across editors with different tabulator widths.

- Structured blocks are indented with four withe spaces. We do not indent the curly braces that
  belong to the structured block. A exception of this rule are the structured blocks of namespaces
  that we do not indent, because this would lead to an unnecessary indention of nearly all our
  source code. The indention of structured blocks greatly improve the readability of the source
  code because it is much easier to spot which lines of code are belonging to the same block.

- Structured blocks are surrounded by curly braces in C-style, i.e. the curly braces **{}** are
  standing below each other in the same column, each in its own row. This rule also applies for one
  statement structred blocks.

Example1

::

   for ( IndexType i = 0; i < n; ++i )
   {
       y[i] += a * x[i];
       if ( doZComputation )
       {
           z[i] = y[i] / b;
       }
   } // for ( IndexType i = 0; i < n; ++i )

Example2

::

   int result = -1;
   switch( a )
   {
       case 0:
       {
           result = 0;
           break;
       }
       case 1:
       {
           result += 1;
           // intended fall through
       }
       case 2:
           // intended fall through
       case 3:
       {
           const int b = 2;
           result += b;
           break;
       }
       default:
       {
           result = -2;
           break;
       }
    }

The arrangement of the curly braces in the same column together improves the ascertainability of
structured blocks.

- Even empty, very short or inline methods should be defined outside the class definition. For
  inline methods and templates this is done in the header file directly below the class definition.
  Besides keeping the interface of the class clean and clear this makes a code rearrangement much
  easier and reduces compile times.

Example

**A.hpp**

::

    class A
    {
    public:
        inline int getI() const;
        int getITimes( int x ) const;
    private:
        int mI;
    };

    int A::getI() const
    {
        return mI;
    }

**A.cpp**

::

    #include "A.hpp"

    int A::getITimes( int x ) const
    {
        return mI * x;
    }

- Whitespaces: We put a blank

   - after each opening bracket
   
   - after each comma
   
   - after each semicolon
   
   - before each closing bracket
   
   - around operators

Naming
^^^^^^

- The upper case letters **I** (i), **O** (o) and lower case letter **l** (L) should not be used alone for an identifier,
  because they can easily be mistaken for an **0** or **1**.

- All Identifiers (class names, function names, variable names, ...) are formatted in CamelCase (Exp.
  **printFunctionName();**). Types are starting with a upper case letter (Exp. **class Matrix;** or **enum
  ExpressionTypes;**) all other identifiers are starting with a lower case letter.

- Makros are all **UPPERCASE_WITH_UNDERSCORE**, should be quite long and should be prefixed with ``LAMA_`` to avoid any
  accidentally replacement by the preprocessor.

- In C source files we do not use CamelCase, because it should be callable by the case insensitiv language Fortran.

Namespace
^^^^^^^^^

- To avoid naming conflicts we use the name space lama in the C++ part and prefix global names with ``lama_`` in the C part.

- The statement **using namespace ...** is not allowed in header files, because it would negate the reason to use name spaces.

- In all none header files **using namespace ...** is consequently used, because it augments the code readability.

Naming conventions
""""""""""""""""""

We should stick to the following naming conventions, because especially when working with sparse matrix formats they greatly enhance the understandability of the source code.

=============== ================================= ===================================================================================================================================================================================
**Name**        **Type**                          **Meaning**
=============== ================================= ===================================================================================================================================================================================
**numRows**     integer                           number of rows in a matrix or a vector
**numColumns**  integer                           number of columns in a matrix
**numValues**   integer                           number of stored elements of a sparse matrix
**ia**          pointer to integer                Index array with row start and end indexes for sparse matrices in compressed row format.
**ja**          pointer to integer                Index array with column indexes for sparse matrices format.
**values**      pointer to floating point number  Floating point array of the values in any matrix of vector format.
**i,j,k**       integer                           Row/column index (in a vector or matrix).
**jj,kk**       integer                           Position in a array, e.g. ja or values. Example **jj** is a position in the array **ja** and points to a column index **j** and the corresponding none zero element in **values**.
=============== ================================= ===================================================================================================================================================================================

To give an example, here the code for a CSR sparse matrix vector multiplication: 

::

   for ( IndexType i = 0; i < numRows; ++i )
   {
       y[i] = 0.0;
       for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
       {
           const IndexType j = ja[jj];
           y[i] += values[jj] * x[j];
       }
   }

Files
^^^^^
- Source files which contain a main method are named like the executable build.

- We use the following file extensions: **.hpp** for C++ header files, **.h** for C header files, **.cuh** for C for
  CUDA header files, **.cpp** for C++ source files, **.c** for C source files and **.cu** for C for CUDA source files.

- All files are named according to their content (class, template, ...) (filename = class name).

- Each file only contains one class or template. This rule should ease the orientation within the project besides that
  smaller files with a single objective lead to less version control conflicts.

- All source code in header files need to use # pragma once. The include header guards avoid violation of the multiple definition rule.
  
- No two files within the project should be only distinguishable through their path or upper and lower case letters.
  This avoids problems with the include header guards and maintains portability.

- Includes

   - Each target directory (lama/, bench/, tests/, ...) has its own system include path, which is set by **include_directories()** in CMake.
    
   - In a target (sub-)directory we use local includes like: **#include "CUDA/CUDADevice.hpp"**
   
   - Avoid using relative paths like **"../Exception/Assert.hpp"**, because in general this indicates a bad design or directory structure.
   
   - **#include <scai/lama.hpp>** from the testing target shouldn't be used. Instead use: **#include <scai/lama/lama.hpp>**.
   
  - Includes of header files that are not part of LAMA should be always done with **#include <file>** to maintain a clear separation between projects.

Example
"""""""
   
directory structure:
   
.. code-block:: bash
   
   src
   src/lama/cppFile.hpp
   src/lama/cppFile.cpp
   src/lama/subdir2/cppFile2.hpp
   src/lama/subdir3/cppFile3.cpp
   src/test/cppFile.hpp
   src/test/cppFile.cpp
   src/test/subdir1/cppFile2.hpp
   
in src/lama/cppfile.cpp:
   
.. code-block:: bash

   #include "cppFile.hpp"
   #include "subdir2/cppFile2.hpp"
   #include <test/cppFile.hpp>
   
in src/lama/subdir3/cppfile3.cpp:
   
.. code-block:: bash

   #include "cppFile.hpp"
   #include "../subdir2/cppFile2.hpp"  //this should be avoided in general
   #include <test/cppFfile.hpp>

in src/test/cppFile.cpp:
   
.. code-block:: bash

   #include "cppFile.hpp"
   #include "cppFile2.hpp"

Directories
^^^^^^^^^^^

- If directories are getting to full sub directories that form logical subgroups should be created.

- Directory names should only consist of lower case ASCII letters to avoid any problems with case
  insensitive file systems (e.g. Windows).

Variables
^^^^^^^^^ 

- We do not use global variables. Global variables make it extremely difficult to spot side effects
  and dependencies and therefore they especially make it hard to parallelize code.

- Variables should be initialized together with their declaration if possible. In general we follow
  the `RAII <http://en.wikipedia.org/wiki/Resource_Acquisition_Is_Initialization>`_ (Resource
  Acquisition Is Initialization) idiom. This helps to avoid memory leaks, especially when exceptions occur.

- Variables should be declared at their first use. This is necessary to follow the RAII idiom and
  restricts the lifetime of the variable and therefore makes the code more maintainable.

- We are `const correct <http://en.wikipedia.org/wiki/Const-correctness>`_, i.e. where it is
  possible to declare a variable as constant we do so. This rule makes it much easier to see
  dependencies and spot side effects.

- Class members are prefixed by **m<VariableName**, e.g. **mIa**. This makes the naming of constructor and getter
  setter methods arguments more easy and helps to determine side effects of methods.

- Class member should be private.

- Static class member are not prefixed with **m**.

- Pointers to accelerator memory (e.g. pointer to the global memory of a CUDA GPU) get the suffix
  **<VariableName>d**, e.g. **mIad**. This avoids confusion between host and device pointers. If it
  is necessary pointers to host memory are have the suffix **<VariableName>h**.

Methods and Functions
^^^^^^^^^^^^^^^^^^^^^

- The argument list of methods and functions are starting with the output parameters. The input (constant) parameters
  follow at the end. Because in an assignment the assigned value is on the left this is intuitive and makes the source
  code more consistent.

Comments
^^^^^^^^

- All comments should conform to the language standard in use, i.e.

  - for comments in C++ source code we use only **//**.
  
  - for comments in C source code we use only **/* */**.
  
- If we have default parameter values or something similar is used we should add a help comment at
  the definition of that function that points to the default parameter. Example **add ( a, b /*=2*/ )**.
  In this case C comments are allowed in C++ source code.

- The closing curly brackets of long structured blocks should be commented with the head of the structured block.

Example:

::

   if ( doLongCalculation )
   {
       // here we have many lines of source code
   } // if ( doLongCalculation )

- All Entities should described with a meaningful doxygen comments at their point of declaration.

  - For doxygen comments C style comments are allowed in C++ source files.
  
  - Example comment for a method

::

	/**
 	* @brief short description
 	* Detailed description
 	* @param([in]|[out]|[in,out]) <parameter-name> parameter description
 	* @throws <class name> exception description
 	* @return description of the return value
 	*/

- Comments within the source code should be inserted if they make the code more readable and easier
  to comprehend. These comments are especially necessary for complicated algorithms, in simpler cases
  describing function and variable names should be preferred. 

Logging and Output
^^^^^^^^^^^^^^^^^^

- We do not use the standard output or standard error streams. Instead we use the logging mechanism
  described in :ref:`logging <scailogging:main-page_logging>`.

- Please keep in mind, that good logging messages are also explaining the source code.

CMake
^^^^^

In all CMake files we stick to the style of the `official CMake documentation`_. All **VARIABLE_NAMES** are written
in upper case and joined by underscore. The **function_names()** are written in lower case letters and are also joined
by underscore. A short introduction to Find modules can be found `here`__.
Modules are named according to the variables they define, e.g. FindSCAI_BLAS and SCAI_BLAS.

.. _official CMake documentation: http://www.cmake.org/cmake/help/documentation.html
__ http://www.itk.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules

.. _tools:

Tools
^^^^^

- Eclipse Code Style:

  - in **tools/eclipse/LAMA-Styleguide.xml** a configuration file for the
    eclipse code formatter is provided. Due to the configurable rules, the configuration file does
    not perfectly fit to this style guide. But it gives a good starting point and therefore it should
    be used.

- astyle:

  - in **tools/lama_format** a shell script which uses the tool astyle to format the source code is provided.

- valgrind:

  - a suppressions file is provided at **tools/valgrind/lama.supp**.

.. _changes:

Changes Version 3.0
===================

Here is a summary of the major changes compared to the old release LAMA 2.1.

 * The build system has been restructured and simplified. There are no more separate build steps 
   for each individual LAMA module. As a build of a module required the dependent modules to be
   installed before, build and installation was mixed up and resulted in many confusions and inconveniences.
   Experience has also shown that individual settings for each module caused more troubles than 
   giving higher flexibility.
 * LAMA now uses full functionality of C++11, therefore the C++ compiler used for LAMA must
   support C++11 features.
 * The classes Matrix and Vector are now template classes with the value type of the elements 
   as template argument. There are still the classes _Matrix and _Vector as common base
   classes, but operations on these base classes are rather limited.
 * In contrary to the previous version, writing code with matrix-vector operations for arbitrary
   types is no longer supported. As a consequence, all solver classes are now also template classes.
 * The class Scalar as common container type for all kind of value types is no longer supported.
   Especially all methods on vector and scalars that return a single value return now a value of 
   the corresponding ValueType and no more a scalar. The class Scalar is still available as it is
   used internally for representation of scalar values in template expressions used to support text book
   syntax, but should not be used any more in user applications itself.
 * Enum types are now all enum classes (C++11 feature). Due to the scoping rules, the qualification of enum values
   is now slightly differrent.
 * The LAMA specific class Thread has been removed and whenever threads are exploited, the corresponding
   classes ``thread``, ``mutex`` and related ones of the C++11 standard are exploited.
 * The support for matrix and vector expressions has improved and is now more functional. Especially unary 
   functions (e.g. sqrt, log, sin, cos ) and binary functions( e.g. min, max, pow) support better text book
   syntax.
 * Type conversions for matrices and vector are still possible but must be called explicitly by a cast function.
 * The LAMA internal types IndexType and PartitionId are now in the namespace scai.
 * The move semantic of C++11 can be used to avoid deep copies of data structures when setting up LAMA
   data structures. This optimizes the code not only for speed but also for need of memory.
 * GPI support (implementation of GASPI standard) has been dropped.
 * MIC support (i.e. support of Xeon Phi Knights Corner in accelerator mode) has been dropped
 * Support of native mode for Xeon Phi Knights Landing has been improved

Typed Vector Classes
--------------------

Methods on matrices and vectors that return a scalar value returned an object of the untyped container
class ``Scalar``. This was due to the fact, that these methods were virtual methods of the base class
``Vector`` and ``Matrix``` that could be used for any type. The inconvienent use of Scalar objects
was also necessary even if the type of the matrix was fixed in the derived class like
``DenseVector<double>``.

.. code-block:: c++

    DenseVector<double> v( 10, 5.0 );
    Scalar alphaS = v[0];
    double v = alphaS.getValue<double>();
    Scalar minS = v.min();
    double minS = minS.getValue<double>();
    Scalar l2NormS = v.l2Norm();
    double l2Norm = l2NormS.getValue<ValueType>;

The classes ``Vector`` and ``Matrix`` are now template classes and so the value type of the elements
is always known.

.. code-block:: c++

    DenseVector<double> v( 10, 5.0 );
    double s = v[0];
    double min = v.min();
    double l2Norm = v.l2Norm();

Polymorphic Code
-----------------

The old LAMA version supported writing polymorphic code without knowing anything about 
the value type at all. 

.. code-block:: c++

    void subroutine( Vector& x, const Matrix& m, const Vector& y )
    {
        VectorPtr z( y.newVector() );
        *z = m * y - x;
        Scalar alpha = z->l2Norm();
        x += alpha * y;
    }

Though this kind of coding was rather general it was not very convenient for 
the user at all. Especially when temporary vectors where required, the use of
pointers was necessary to create vectors of the same value type. Even if the use
of different value types in expressions was allowed, it was nearly impossible to 
identify where this was the case and what kind of support or overhead was implied.

With the latest version, this code has to be rewritten as template code to keep
the generality.

.. code-block:: c++

    template<typename ValueType>
    void subroutine( Vector<ValueType>& x, const Matrix<ValueType>& m, const Vector<ValueType>& y )
    {
        DenseVector<ValueType> z;
        z = m * y - x;
        ValueType alpha = z.l2Norm();
        x += alpha * z;
    }

There is no more need to use a pointer variable for the temporary vector as now the type
of the vector is explicitly given. Furthermore the use of the class Scalar is no
more necessary and the used ValueType is now always explicit.

Matrix-Vector Expressions
-------------------------

The support of text-book syntax has been improved.

.. code-block:: c++

   DenseVector<double> x;
   DenseVector<double> y;

   x.setRandom( 1000, 1 );    // 1000 random value between 0 and 1
   x = sin ( x );             // was x.sin() before
   x = cos ( x );             // was x.cos() before
   y = exp ( x );             // same as y = x; y.exp() before

   y = 1 / y;                 // y.invert() before;
   y = x / y;                 // was not supported before


While in the previous release mixing of different value types was nearly
always possible, this is now rather restricted.

.. code-block:: c++

     DenseVector<float> x( 100, 1.0f );
     DenseVector<double> y( 100, 2.0 );
     DenseVector<float> z;
     z = x + y;       // you cannot mix different types in expressions

     z = cast<float>( y );   // conversion of double to float is supported.
     z = x + z;              // this is now an expression with float operands

Beside some exceptions, in the old version the use of mixed value types was not well
supported and as a result temporary vectors where created. 

.. code-block:: c++

   z = x + y;   // requires temporary vector for y

   DenseVector<float> tmpY = y;
   z = x + tmpY;

The solution with using z itself as temporary version for the conversion of y is much more
efficient than using a new temporary vector. 

Type Conversions
----------------

Format and type conversions of matrices and vector were supported well by previous LAMA releases.
This functionality is still available but must be called explicitly. 

.. code-block:: c++

   DenseVector<float> x( "input.txt" );
   DenseVector<double> z( x );          // okay, type conversion in copy constructor
   SparseVector<double> y( x );         // okay, format and type conversion in copy constructor
  
   z = x;                    // is now illegal, was supported in previous versions
   z = cast<dense>( x );     // that is now the right way for type conversion
   y = cast<double>( x );    // type conversion can also involve format conversion
   z = y;                    // format conversions are still done implicitly.

Type conversions in expressions are no more supported.

.. code-block:: c++

   DenseVector<float> fD;
   DenseVector<double> dD;
   SparseVector<float> fS;
   SparseVector<double> dS;
   ...
   fD = 2 * fD + fS;                   // okay, all vectors are float
   dD = 2 * fD + dS;                   // ERROR, cannot mix float and double
   dD = 2 * cast<double>( fd ) + dS;   // ERROR, cast within expression, might imply use of temporary vectors
   dD = cast<double>( fd );            // okay, conversion is done in-place
   dD = 2 * dD + dS;                   // fine, reuse of dD avoids a temporary

Enum Classes
------------

In the old release plain enums where exploited, but in an own name space as the following example shows:

.. code-block:: c++

    namespace scalar {
  
        enum ScalarType = { REAL, FLOAT, ... };
    }

The new release now uses enum classes and so the own namespace is no more required. 

    enum class ScalarType = { REAL, FLOAT, ... };

Unfortunately this causes some renaming, but it avoids confusion between the name of the namespace and the
name of the enum type.

.. code-block:: c++

     scalar::ScalarType s = scalar::REAL;   // old code
     ScalarType s = ScalarType::REAL;       // new code

Simliar other examples are:

.. code-block:: c++

    common::binary::BinaryOp op = common::binary::COPY;
    common::BinaryOp op = common::BinaryOp::COPY;  

Shared and Unique Pointers
--------------------------

In the previous release LAMA provided own classes for shared and unique pointers.
Actually these classes were wrappers for the std classes shared_ptr and unique_ptr
or for the corresponding Boost classes (for compilers without C++11 support).

As LAMA now relies on C++11 support, these wrapper classes became redundant.
They have been removed completely as many user applications are using these 
pointer classes already for themselves and the functionality of these wrapper
classes was slightly limited (e.g. unique_ptr could not be used in C++ container
classes like vector).

.. code-block:: c++

    #include <scai/common/unique_ptr.hpp>
    #include <scai/common/shared_ptr.hpp>

    using namespace scai;

    common::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    common::shared_ptr<lama::Vector> vzzPtr( vZ.newVector() );
    common::scoped_array<double> mG( new double[10] );

The changes required for the new LAMA version are rather straightforward.

.. code-block:: c++

    #include <memory>

    std::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    std::shared_ptr<lama::Vector> vzzPtr( vZ.newVector() );
    std::unique_ptr<double[]> mG( new double[10] );

The use of  the pointer variables itself does not require any changes.

Move Semantic
-------------

Here is a typical LAMA code of how to set up a CSR sparse matrix with raw data.

.. code-block:: c++

    // copy the raw data into heterogeneous array
    HArray<IndexType> csrIA( numRows + 1, rawIA );
    HArray<IndexType> csrJA( numValues, rawJA );
    HArray<ValueType> csrValues( numValues, rawValues );
    // build a CSR storage, copies the input arrays
    CSRStorage<ValueType> csrStorage ( numRows, numColumns, csrIA, csrJA, csrValues );
    // build a CSR matrix, copies the CSR storage
    CSRSparseMatrix<ValueType> csrMatrix( csrStorage );

One copy, here the first one, is mandatory as otherwise data cannot be modified and managed on its own.
The two other copies are not really needed. Actually we only want to move the allocated data of the containers.
By using the move semantics of C++11, it is possible to avoid these copy steps as follows:

.. code-block:: c++

    CSRStorage<ValueType> csrStorage ( numRows, numColumns, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );
    CSRSparseMatrix<ValueType> csrMatrix( std::move( csrStorage ) );

Please note that the move operations leave the heterogeneous arrays ``csrIA``, ``csrJA``, and ``csrValues``
as well as the CSR storage ``csrStorage`` as empty containers.

By using the move operations the code becomes faster and it might be helpful to avoid running out of memory.

AssemblyAccess
--------------

The class SparseAssemblyStorage supported in former LAMA releases is no more available.

Instead of this you shoud use the class	MatrixAssemblyAccess that has a similiar functionality
but works also well for distributed matrix data.

The class VectorAssemblyAccess has a similiar functionaly fot distributed vectors.

Tutorial/Lecture
-----------------

The tutorial and lecture in the user guide have been completely revised.

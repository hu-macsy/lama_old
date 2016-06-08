Class LArray
============

The array class **LArray** stands for a local array and is derived from the class *HArray*. Regarding
the member variables it is exactly the same, but it offers many addtional operators.

.. code-block:: c++

    template<typename T>
    class LArray : HArray<T>
    { 
    }

Added Operations
----------------
 
- constructors
 
  - copy constructors for LArray and _HArray
   
  - initialization with number of entries and a normal array
 
- operators ( usage with single value or another _HArray )
 
  - operator=
   
  - operator*=
 
  - operator/=
 
  - operator+=
 
  - operator-=
 
  - operator[]
 
- reduction
 
  - min
 
  - max
 
  - sum
 
- transform
 
  - invert
 
  - conjuagte
  
- calculation
  
  - dot-product
 
- norm
 
  - L1-norm
 
  - L2-norm
 
  - Max-norm
 
  - Max-Diff-norm
  
Usage
-----

For the implementation of the operators kernel functions implemented on different devices
are used.

.. code-block:: c++

    LArray<double> A;
    input( A );
    LArray<float> B( A );   // implicit conversion
    ...

    A = A * B;   // componentwise multiplication  
    
    A *= 2.0; // scaling
    
    B -= 13; // componentwise subtraction
    
    double d = A.dotProduct( B ); // vector multiplication
    
    float mi = B.min(); // minimum
    
    double ma = A.max(); // maximum
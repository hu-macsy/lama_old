KernelTrait
===========

A ``KernelTrait`` is a struct that has an entry for the name of the kernel routine and a type 
definition for its signature.


.. code-block:: c++

   struct function_name
   {
       typedef res_type ( *FuncType ) ( type1, ..., typen );
       static const char* getId() { return "function_name"; }
   };


Template arguments can be used to define traits for different value types.


.. code-block:: c++

   template <typename T1, typename T2, typename T3>
   struct template_function_name
   {
       typedef res_type ( *FuncType ) ( T1*, const T2*, T3 );
       static const char* getId() { return "template_function_name"; }
   };
   

A ``KernelTrait`` is a struct that is introduced for the following reasons:

- it avoids misspelling of function names

  - different strings used for registration and access will not be detected at compile time, but different struct identifiers.

- The same is true for the correct signature, i.e. functions are never registered and used with different signatures.

- They are used for Doxygen documentation of the function behavior


Multiple ``KernelTrait`` might be grouped for kernel routines used in a certain module. The 
name of the group should appear as prefix in the id of the kernel routines.

.. code-block:: c++

  struct UtilKernelTrait
  {
      template <typename ValueType>
      struct maxval
      {   /** @brief Find maximal value of n contiguously stored values.
           *  @param[in] array is an array of values
           *  @param[in] n is the size of array
           *  @return maximum of all values in array                      
          */
          typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );

          static const char* getId() { return "Utils.maxval"; }
      };
  
      // Other traits of Utils routines
  
      ...

  }; // struct UtilKernelTrait
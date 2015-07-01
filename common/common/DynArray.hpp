#pragma once

/** \file DynArray.hpp    

    \brief Definition of a template classes for dynamic arrays
           (up to 5 dimensions).
*/

#include <sstream>
#include <stdexcept>

namespace common {

  /**************************************************************************
  *                                                                         *
  *  DynArray1  - one dimensional dynamic arrays                            *
  *                                                                         *
  **************************************************************************/

  /** This class maintains dynamic one-dimensional arrays of any type.

      A one-dimensional array will be allocated in a contiguous memory 
      array of size n. 

      * Initialization of array elements is not supported
      * copy constructor is disabled

      \code
      void calc(double* m);
      DynArray1<double> mass(nbody, "Module::mass");
      for (int i = 0; i < mass.size(); i++) {
        mass[i] = 0.0;
      }
      calc(mass);   // type conversion to double* supported
      mass.free();  // usually done only implicitly by destructor
      \endcode
  */

  template <typename T> 
  class DynArray1 {

    public:

      /** Default constructor creates a NULL array. */

      DynArray1()
      {
        data = NULL;
        n = 0;
      }

      /** Construct a dynamic array with given size.

          \param n is the size of the array
          \param msg is a string with some info about the array
      */

      DynArray1(size_dyn n, const char* msg)
      {
        data = NULL;  // important before allocate, avoids free of illegal ptr
        allocate(n, msg);
      }

      /** Destructor will free the allocated memory. */

      ~DynArray1()
      {
        if (data) {
          free();
        }
      }

      /** Assignment operator for initializations. */

      DynArray1<T>& operator = (const T val)
      {
        for (size_dyn i = 0; i < n; i++) data[i] = val;
        return *this;
      }

      /** Getter for the size of the array. */

      size_dyn size() const {
        return n;
      }

      /** Allocate the array with a given size. Works also as reallocate. 

          \param n is the size of the array
          \param msg is a string with some info about the array
      */

      void allocate(size_dyn n, const char* name) {

        if (data != NULL) Memory::sfree(data);
        data = (T*) Memory::smalloc(n * sizeof(T), name);
        this->n = n;
      }

      /** Explicitly free the memory. Will be called in any
          case when the object is destroyed.
      */

      void free()
      {
        Memory::sfree(data);
        data = NULL;          // Important to nullify
        n = 0;
      }

      /** Extend array in first dimension, existing values are not saved
       *
       * \param n is the new size
       * \param msg is used as log/error message
       *
       * Extension of the array will only be done if the current size is not
       * large enough.
       *
       * */

      void extend(size_dyn n, const char* msg)
      {
        if (n <= this->n) return;  // sufficient
        allocate(n, msg);
      }

      void extend(size_dyn n, size_dyn nmax, const char* name)
      {
        if (n <= this->n) return;
        allocate(nmax, name);
      }

      /** extend array in first dimension and keep old values */

      void grow(size_dyn n, const char* name)
      {
        if (n <= this->n) return;
        data = (T*) Memory::srealloc(data, n * sizeof(T), name);
        this->n = n;
      }

      const T& operator[] (size_dyn i) const {
        return  data[i];
      }
 
      T& operator[] (size_dyn i) {
        return  data[i];
      }

      /** Get reference to an element but with index check. */

      T& at(size_dyn index)
      {
        // Note: comparison < 0 not needed as it is an unsigned int

        if (index >= n) {
          throw std::out_of_range("DynArray1::at");
        }
        return data[index];
      }

      /** Get const reference to an element but with index check. */

      const T& at(size_dyn index) const
      {
        if (index >= n) {
          throw std::out_of_range("DynArray1::at");
        }
        return data[index];
      }

      /** Query operator if array has been allocated.

          \return true only if memory has been allocated

          Note: returns also false, if array has been allocated with
          at least one zero dimension.
      */

      operator bool() const {
        return data != NULL;
      }

      /** memory usage of the dynamic array

          \return number of allocated bytes
      */

      size_dyn memoryUsage() const
      {
        return n * sizeof(T);
      }

      /** Conversion operator to get pointer to the data.

          \code
          void jacobi(double* matrix);
          ...
          DynArray<int> x(n1, n2);
          jacobi(x);
          \endcode
      */

      operator T*() { return data; }

      /** Conversion operator to get const pointer to the data.

          \code
          double sumit(const double* matrix);
          ...
          DynArray<int> x(n1, n2);
          double a = sumit(x);
          \endcode
      */

      operator const T*() const { return data; }

      /** Switching contents of two dynamic arrays. */

      void swap(DynArray1<T>& other) 
      {
        T* swapPtr = other.data;
        other.data = data;
        data = swapPtr;

        size_dyn swapN = other.n;
        other.n = n;
        n = swapN;
      }

    private:

      size_dyn n;      //!< allocated size
      T* data;         //!< pointer to the first element
   
      // disable the copy constructor

      DynArray1(const DynArray1&) {}

  };

}


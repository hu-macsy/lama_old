.. _Math:

Math
====

Math is a wrapper for functions from cmath. Due to the compatibility of cmath to C it is not the case that mathematical operations are overloaded
for all arithmetic types. E.g. there are different routines like sqrtf, sqrt and sqrtl
for the square root function. The Math wrapper enables you
just to write Math::sqrt. This is very helpful to write templated code. Besides sqrt Math
provides wrapping functions for abs and conj. Currently int, long, long, float, double, long double,
Complex<float>, Complex<double> and Complex<long double> are supported.  

The following template function returns the L2-norm for an array that contains
values of a certain type. If you call this code with ValueType=float it will lead to compiler warning, 
because the sqrt-function expects an argument of type double. 

.. code-block:: c++

  template<typename ValueType>
  ValueType l2norm( const ValueType vals[], const int n )
  {
      ValueType sum = 0;
      for ( int i = 0; i < n; ++i )
      {
          sum += vals[i] * vals[i];
      }
      return sqrt( sum );
  }


For this and some other situations, a struct Math is used that provides 
type-specific functions or values that can be used. In the example above, the 
code becomes:

.. code-block:: c++

  template<typename ValueType>
  ValueType l2norm( const ValueType vals[], const int n )
  {
      ValueType sum = 0;
      for ( int i = 0; i < n; ++i )
      {
          sum += vals[i] * vals[i];
      }
      return Math::sqrt( sum );
  }

Beside the mathematical routines, the Math class provides also a common facitlity to generate
random numbers.

* ``Math::random<ValueType>( unsigned nb )`` return a value between 0 and nb for any value type
* ``Math::randomBool( const float trueRatio )`` generates boolean value where the value true has of a probability of trueRatio.
* ``Math::srandom( int )`` sets the seed for the randon number generator, allows to reproduce results or to force different results on different processors with different seeds.


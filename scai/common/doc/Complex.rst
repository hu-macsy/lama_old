.. _Complex:

Complex
=======

The common library provides an own implementation for complex numbers instead
of using the std::complex type. The template class Complex can be instantiated
for the types float, double or long double. Here is a short example:

.. code-block:: c++

  #include <scai/common/Complex.hpp>

  using scai::common::Complex;
 
  Complex<float> a( 1.0f, 1.0f );
  Complex<double> b( 2.0, -1.0 );
  std::cout << "a = " << a << ", b = " << b << std::endl;


Advantages:

 * output operator just prints the two values re and im
 * provides conversion operators needed for I/O
 * Implements fabs, operator<, and operator>
 * can also be used in CUDA.

Typedefs are used to have one single name for the different instantiations.

.. code-block:: c++

   typedef scai::common::Complex<float> ComplexFloat;
   typedef scai::common::Complex<double> ComplexDouble;
   typedef scai::common::Complex<long double> ComplexLongDouble;

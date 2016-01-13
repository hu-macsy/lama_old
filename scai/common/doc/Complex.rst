Complex
-------

The common library provides an own implementation for complex numbers instead
of using the std::complex type.

.. code-block:: c++

    namespace scai
    {
        namespace common
        {
            template<typename ValueType>
            class Complex
            {
                ...
            }
        }
    }

Advantages:

 * provides conversion operators needed for I/O
 * Implements fabs, operator<, and operator>
 * can also be used in CUDA kernel implementations.

Typedefs are used to have one single name for the different instantiations.

.. code-block:: c++

   typedef scai::common::Complex<float> ComplexFloat;
   typedef scai::common::Complex<double> ComplexDouble;
   typedef scai::common::Complex<long double> ComplexLongDouble;

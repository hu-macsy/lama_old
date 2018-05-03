Free Functions for HArrays
===========================

Some functions are provided to construct objects of the class HArray.

.. code-block:: c++

    // create an array of size 100 filled with 0
    auto zeroArray = utilskernel::fillHArray<double>( 100, 0.0, ctx );

    // create an array of size 100 filled with random values between 0 and high 
    auto randomArray = utilskernel::randomHArray<double>( 100, high, ctx );

    // create an array of size 100 sparsely filled, i.e 20% are filled with values between 0 and high
    // while the other values are set to 0 ( where the default can be any value )
    auto sparseArray = utilskernel::sparseRandomHArray<double>( 100, 0, 0.2f, high, ctx );

    // convert a double array to a float array
    auto sparseFArray = utilskernel::convertHArray<float>( sparseArray );


First LAMA Example Program
==========================

The following C++ program shows a very simple example program of how to use LAMA.

::

    // Always include this file at the beginning

    #include <lama.hpp>

    #include <lama/DenseVector.hpp>
    #include <lama/Scalar.hpp>

    #include <iostream>

    using namespace lama;

    int main()
    {
        IndexType size = 8;
        DenseVector<double> v( size, 1.1 );
        Scalar s = v.l1Norm();
        std::cout << "L1 norm of v = " << s.getValue<double>() << std::endl;
        return 0;
    }

:download:`Download source file <../../cpp_source/tutorial/simple.cpp>`

The include file lama.hpp contains some definitions how far assertions, logging and tracing statements
are compiled into your code. The definitions will be the same as used for the installation.

Usually you have to include the class definition file for each LAMA class that you are
using. As we use objects of class DenseVector and Scalar, we have to include the corresponding files.


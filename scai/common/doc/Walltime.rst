.. _Walltime:

Walltime
========

The class Walltime provides a portable routine to get the current walltime in seconds.
The following example shows a typical use how to measure the walltime spent in a certain
piece of code:

.. code-block:: c++

  #include <scai/common/Walltime.hpp>

  double time = scai::common::Walltime::get();

  // code to measure walltime

  ...

  time = ( scai::common::Walltime::get() - time );

  cout << "Code needed " << time << " seconds" << endl;

Furthermore, the class provides a portable sleep routine that makes the calling thread sleeping for a 
certain number of milliseconds.

.. code-block:: c++

  common::Walltime::sleep( 500 );  // sleep 500 milliseconds

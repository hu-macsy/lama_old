Timing
------

.. literalinclude:: ../Walltime.hpp
   :language: c++
   :lines: 53-78

.. code-block:: c++

  double time = scai::common::Walltime::get();

  // code to measure walltime

  ...

  time = ( scai::common::Walltime::get() - time );

  cout << "Code needed " << time << " seconds" << endl;



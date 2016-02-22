Timing
------

.. code-block:: c++

  typedef uint64_t INTEGER_8;

  class COMMON_DLL_IMPORTEXPORT Walltime
  {
  public:

      /**  @return current walltime in ticks */
      static INTEGER_8 timestamp();
 
      /** Number of ticks per second, so timestamp() / timerate() gives time in seconds. */
      static INTEGER_8 timerate();
  
      /** @return current walltime in seconds */
      static double get();
  };

.. code-block:: c++

  double time = scai::common::Walltime::get();

  // code to measure walltime

  ...

  time = ( scai::common::Walltime::get() - time );

  cout << "Code needed " << time << " seconds" << endl;



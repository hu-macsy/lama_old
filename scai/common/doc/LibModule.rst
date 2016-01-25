LibModule
---------

Many functionality of LAMA can be 
either linked directly with the LAMA library
Or can be loaded dynamically after runtime.

Dynamic extension is supported by:

 * Factory (add new Context, Memory, ...)
 * KernelRegistry (add new kernel routines)

Advantages of using Library Modules loaded at runtime:

 * They are really loaded (for linking –Wl,no-as-needed is required)
 * Order of modules can be used to overwrite entries in registries

Be careful: unload/free of module not mandatory

.. code-block:: c++

  class LibModule
  {
  public:

      /** Data type definition for library handle, might be OS specifici. */

      typedef void* LibHandle;

      /** Load a library with its full name (might be absolute or relative) */

       static LibHandle loadLib( const char* filename );

       /** Unload a library  */

       static void freeLib( LibHandle handle );

       /** This routine reads a directory and tries to load all library modules in it */

      static void loadLibsInDir( const char* dir );

      /** multiple directory/libaries by patterns <pattern1>:<pattern2>:<pattern3>  */

      static void loadLibsByPath( const char* path );
  }

Example of a dynamic module:

.. code-block:: c++

  class DynRoutine  : public scai::common::Factory<std::string, DynRoutine*> 
  {
  public:

      /** This routine must be provided by all derived classes. */

      virtual void doIt() = 0;
  };

  class Function1 : public DynRoutine, public DynRoutine::Register<Function1>
  {
  public:
  
      /** This is the implementation for the virtual function of base class DynRoutine */

      void doIt() { ... }

      /** Value for which this class registers in the factory. */

      static inline std::string createValue() { return "Function1";  }

      /** Method that creates objects of type Function1 that will be used for registration. */

      static DynRoutine* create() { return new Function1(); }
  };

Example of using a dynamic module:

.. code-block:: c++

  void runIt()
  {
      if ( !DynRoutine::canCreate( "Function1" ) )
      {
          std::cout << "Function1 not available" << std::endl;
          return;
      }

      std::cout << "Function1 can be created." << std::endl;

      unique_ptr<DynRoutine> d ( DynRoutine::create( "Function1" ) );

      d->doIt();
  }

.. code-block:: c++

   runIt();  // Function1 not available

   LibModule::LibHandle handle = LibModule::loadLib( lib );

   runIt();  // Function1 can be created

   LibModule::freeLib( handle );  // unregister now supported by Factory

   runIt(); // Function1 not available


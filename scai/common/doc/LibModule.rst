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

..
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
	  
.. literalinclude:: ../LibModule.hpp
   :language: c++
   :lines: 47-78


Example of a dynamic module:

.. literalinclude:: ../examples/DynRoutine.hpp
   :language: c++
   :lines: 22-
   
.. literalinclude:: ../examples/Module.cpp
   :language: c++
   :lines: 10-34

Example of using a dynamic module:
	  
.. literalinclude:: ../examples/UseModule.cpp
   :language: c++
   :lines: 13-28	  

.. code-block:: c++

   runIt();  // Function1 not available

   LibModule::LibHandle handle = LibModule::loadLib( lib );

   runIt();  // Function1 can be created

   LibModule::freeLib( handle );  // unregister now supported by Factory

   runIt(); // Function1 not available


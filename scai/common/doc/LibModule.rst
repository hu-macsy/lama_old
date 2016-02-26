LibModule
=========

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


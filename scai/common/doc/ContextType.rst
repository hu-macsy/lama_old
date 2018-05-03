.. _ContextType:

ContextType
===========

A context specifies a device like CPU or GPU accelerator where memory can be
allocated and where operations on it might be executed.
The enumeration type for context might be used for registration and 
searching of objects belonging to a certain context.

===============        =============================
Name                   Meaning
===============        =============================
``Host``               context for cpu + main memory
``CUDA``               CUDA GPU device
``OpenCL``             OpenCL GPU device, currently not supported
``UserContext``        can be used for a new derived Context class
``MaxContext``         dummy value, used for dimension of ContextType arrays
===============        =============================

Note: Extending LAMA with support of a new context implies to add an item to this enum type.


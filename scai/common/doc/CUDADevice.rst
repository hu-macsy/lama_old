CUDADevice
==========

The class CUDADevice takes core for the access of a CUDA device. It checks for
the availability and creates a CUDA context structure that is required for
further API calls of CUDA routines. The destructor of the object destroys the
context.

When the first CUDADevice object is created, CUDA will also be initialized.

.. literalinclude:: ../examples/cuda/CUDADeviceExample.cpp

In the above example a CUDA device is created where the number of the device
can have been set via the environment variable ``SCAI_DEVICE``.

 * the getter method getCUcontext returns the CUDA context structure
 * the getter method getCUdevice returns the CUDA device structure

Before any API routine of the CUDA driver or CUDA runtime is called 
that requires a CUDA context argument, it should be made
sure that the corresponding CUDA context is pushed before. 

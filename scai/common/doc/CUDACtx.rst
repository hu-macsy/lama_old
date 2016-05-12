.. _CUDACtx:

CUDACtx
=======

The class CUDACtx takes care for the access of a CUDA device. It checks for
the availability and creates a CUDA context structure that is required for
further API calls of CUDA routines. The destructor of the object destroys the
context.

When the first CUDACtx object is created, CUDA will also be initialized.

.. literalinclude:: ../examples/cuda/DemoCUDACtx.cpp

In the above example a CUDA context is created where the number of the device
can be set via the environment variable ``SCAI_DEVICE``.

 * the getter method getCUcontext returns the CUDA context structure
 * the getter method getCUdevice returns the CUDA device structure

Before any API routine of the CUDA driver or CUDA runtime is called 
that requires a CUDA context being enabled, it should be made
sure that the corresponding CUDA context is pushed before. 

It is possible to use multiple CUDA devices in one application each one with
its own context. Even the same CUDA device might be used logically twice in one application
by creating two contexts.

.. code-block:: c++

    CUDACtx context1( 0 );
    CUDACtx context2( 1 );
    CUDACtx context3( 0 );   

So even if context1 and context3 use the same physical device 0, they behave
like two different devices. Memory allocated for context1 cannot be used
directly on context3 and vice versa. The support of multiple devices is
also one reason why a device must be accessed explicitly by pushing the
corresponding context. Otherwise it is undefined on which device e.g.
a memory allocation or a kernel launch is done. 

Note: LAMA itself will always use only one context on one device.


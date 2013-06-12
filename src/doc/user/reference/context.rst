Contexts
========

A *Context* is concept used by LAMA to define where to execute a calculation and therefor where to store the data.
The default context is always the host (CPU). Beneath the host context there exists two other contexts in LAMA: CUDA and
OpenCL (in progress), which are located on the GPU (for CUDA) and other OpenCL supported accelertors (for OpenCL).

Contexts are registered by a factory

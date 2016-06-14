.. _MIC:

Using Intel MIC in LAMA
^^^^^^^^^^^^^^^^^^^^^^^

MIC stands for many integrated core. The Intel Xeon Phi is a device that has this architecture. 

LAMA suports the Intel MIC as a device on which data can be allocated and on which compute kernels can
be executed. Data transfer between Host CPU memory and Intel MIC memory is also supported. Therefore
it can be used in a very similiar way like a NVIDIA GPU. 

The current release uses for the implemenation of kernels and other operations on the Intel MIC the offload programming model.
Therefore an Intel C++ compiler is required that supports the corresponding offload pragmas. 
These compilers implement also the OpenMP API that is used to benefit of the many cores on one device.
Furthermore, the Intel MKL 
version for the Intel MIC can be exploited for the implemenation of many kernel routines. This is not  mandatory but highly
recommended to achieve better performance if these routines are cruicial for the performance of your application.

The support in LAMA is an optional feature. To enable it you have to pass the parameters ``-DUSE_MIC=1`` to cmake. 
The configuration will fail if the chosen C++ compiler does not support the offload model.

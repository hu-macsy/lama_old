Using Intel MIC in LAMA
^^^^^^^^^^^^^^^^^^^^^^^

The MIC stands for many integrated core. The Intel Xeon Phi is device which has this architecture. 

To utilize it a corresponding Intel compiler and the MKL is needed. Currently the mic backend is parallized by using OpenMP.  

The support in LAMA is an optional feature. To enable it you have to pass the parameters -DUSE_MIC=1 to cmake. 

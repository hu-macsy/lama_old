Examples for calls
==================
matrixGenerator.exe example 3 27 100 100 100

creates a 3-dimensional 27-point poisson star with dimension 100 in every direction
(supported poisson stars: 1D3P 2D5P 2D9P 3D7P 3D27P)

1 CPU:   lamaSolver.exe example [ --SCAI_CONTEXT=Host|CUDA ]
1 GPU:   lamaSolver.exe example CUDA [ ELL ]
1 GPU:   lamaSolver.exe example CUDA JDS
2 CPU:   mpirun -np 2 lamaSolver.exe example 
2 GPU:   mpirun -np 2 lamaSolver.exe example CUDA JDS

mpirun -bind-to-socket --map-by socket -np 2 lamaSolver.exe ...

Further arguments:

   Wxxx gives the process a weight,

        mpirun -np 2 lamaSolver.exe ... W1,W2

        first process gets approximately half number of rows as second process

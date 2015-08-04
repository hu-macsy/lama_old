Examples for calls
==================
matrix_generator.exe example 3 27 100 100 100

creates a 3-dimensional 27-point poisson star with dimension 100 in every direction
(supported poisson stars: 1D3P 2D5P 2D9P 3D7P 3D27P)

1 CPU:   cg_solver.exe example [ Host CSR }
1 GPU:   cg_solver.exe example CUDA [ ELL ]
1 GPU:   cg_solver.exe example CUDA JDS
2 CPU:   mpirun -np 2 cg_solver.exe example 
2 GPU:   mpirun -np 2 cg_solver.exe example CUDA JDS

mpirun -bind-to-socket -bysocket -np 2 cg_solver.exe ...

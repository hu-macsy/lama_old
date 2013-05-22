Examples for calls
==================

matrix_generator.exe example 1 3 1000000
matrix_generator.exe example 2 5 1000 1000
matrix_generator.exe example 2 9 1000 1000
matrix_generator.exe example 3 7 100 100 100
matrix_generator.exe example 3 27 100 100 100

1 CPU:   cg_solver.exe example [ Host CSR }
1 GPU:   cg_solver.exe example CUDA [ ELL ]
1 GPU:   cg_solver.exe example CUDA JDS
2 CPU:   mpirun -np 2 cg_solver.exe example 
2 GPU:   mpirun -np 2 cg_solver.exe example CUDA JDS

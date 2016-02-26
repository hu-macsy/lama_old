Adding a new Sparse Matrix Storage
----------------------------------

- Add MatrixStorage Class (See CSRStorage as an Example)

- Add Matrix template by deriving from the template MatrixTemplate (see CSRSparseMatrix as an example)

  - Add a Test PMatrixTest to test serial and parallel construction
  
- Add OpenMP Interface for the new SparseMatrix format ( see OpenMPCSRInterface as an example)

- Add Backend Interfaces (CUDA, OpenCL, ...) for the new SparseMatrix format ( see CUDACSRInterface for CUDA
  as an example )

- Test Matrix Vector Multiplication by instantiating matrixTimesVectorTestImpl for the new Matrix template in
  PVectorTest

- optional: Add a Benchmark for Matrix Vector Multiplication by instantiating LAMAMVBenchmark for the new Matrix template

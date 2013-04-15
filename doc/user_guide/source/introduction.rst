************
Introduction
************

This documentation will give an overview of the features of the \Library of
\Accelerated \Math \Applications (\L\A\M\A) and will describe the general usage
of these features.

You will find concise :ref:`installation_index` instruction, an easy beginners guide with serval short describtions of
our basic data structures, their behaviour as well as simple examples on how to use them (:ref:`tutorial_index`). An
advanced :ref:`lecture_index` of our hands-on session follows before you get the chance to have a look at the
:ref:`reference_index`. We close with informations about writing own tests (see :ref:`testing_index`) and 
:ref:`benchmarks_index` and on special hints :ref:`developer_index`. 

About LAMA
==========

.. image:: _images/LAMA.png

LAMA is an easy to use open source \Basic \Linear \Algebra \Subprogram (:ref:`\B\L\A\S <blas_explanation>`)) Library with special focus on large
sparse matrices. Its core is written in C++, so you have the comfort of writting your algorithms in a text-book syntax
because of expressions templates. Within the given backends there are optimized matrix-vector kernels for various sparse
matrix formts organized, so you do not have to care about hardware specific programming (e.g. CUDA). We also care about
communicating distributed data between processes. On top of this we prepare a set of basic linear solvers as a jacobi
or CG method.

Our goals
=========

- easy to use text-book syntax for linear algebra operations (y = A * x + b)
- hidden complexity of hardware specific programming and communication needs 
- easy integration of different heterogeneous hardware components, especially GPUs
- high parallelism ( in core concurrency and multinode parallelism)
- easy extensibility (e.g. for sparse matrix formats, compute backends, communication models, ...)

Given features
==============

The supported features of the actual release are listed below:

+------------------------------+----------------------------------------------------+
| - matrix formats             | - distributions                                    |
|                              |                                                    |
|   - Dense                    |   - blocked                                        |
|                              |                                                    |
|   - Sparse                   |   - cyclic                                         |
|                              |                                                    |
|     - CSR                    |   - blockcyclic (general block)                    |
|                              |                                                    |
|     - COO                    |   - general                                        |
|                              |                                                    |
|     - ELL                    | - solvers                                          |
|                              |                                                    |
|     - JDS                    |   - direct solvers                                 |
|                              |                                                    |
|     - DIA                    |     - inverse solver                               |
|                              |                                                    |
|  - backends                  |   - iterative solvers                              |
|                              |                                                    |
|    - CPU (OpenMP optimized)  |     - jacobi method                                |
|                              |                                                    |
|    - GPU                     |     - conjugated gradiant (CG) method              |
|                              |                                                    |
|      - CUDA                  |     - generalized minimal residual (GMRES) method  |
|                              |                                                    |
|                              |     - simple algebraic multigrid (SAMG) method     |
|                              |                                                    |
|                              |     - successive over-relaxation (SOR) method      |
+------------------------------+----------------------------------------------------+
    
Work in process that you can find in our feature branches and that will come with next releases are:
 
- OpenCL-Backend

- C-Interface

- PGAS support for communication (till now we support MPI)

- ease to use configuration of solver through a DSL

- sparse matrix ordering and partitioning through METIS

- connectivity to OpenFOAM

There are also a couple of features that are planned:

- consideration of GPU direct

- OpenACC backend

- MatLab interface

- support of structured matrices

- mixed precision

- complex data type

Additionally to our open source implementation we also support a commercial version of the AMG method.
If you are interested in this, please `contact`_ us. 

.. _`contact`: mailto:lama@scai.fraunhofer.de

Use Cases
=========

So, when it's time to use LAMA?

LAMA is the right decision for you, if you are doing linear algebra on sparse matrices and you want to reach the full
performance of your (parallel) machine without taking care on the kernel code on your own. With LAMA it's easy for you
to write code, that is executable on different heterogeneous machines. 

Possible use cases for LAMA can be found in solving partial elliptic differential equation (e.g. for ground water flow
and oil reservation simulation), image filtering (blur, gauss filtering) and many more. 

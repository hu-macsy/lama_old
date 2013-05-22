************
Introduction
************

This documentation will give an overview of the features of the \Library of \Accelerated \Math \Applications
(\L\A\M\A) and will describe the general usage of these features.

You will find concise :doc:`installation/index` instruction, an easy beginners guide with several short descriptions of
our basic data structures, their behavior as well as simple examples on how to use them (:doc:`tutorial/index`). An
advanced :doc:`lecture/index` of our hands-on session follows before you get the chance to have a look at the
:doc:`reference/index`. We close with informations about writing own tests (see :doc:`testing/index`) and 
:doc:`benchmarks/index` and on special hints :doc:`developer/index`. 

About LAMA
==========

.. image:: _images/LAMA.png
   :align: center
   :alt: LAMA Design

LAMA is an easy to use open source \Basic \Linear \Algebra \Subprogram (:doc:`\B\L\A\S <blas_explanation>`) Library with
special focus on large sparse matrices. Its core is written in C++, so you have the comfort of writing your algorithms
in a text-book syntax because of expressions templates. Within the given backends there are optimized matrix-vector
kernels for various sparse matrix formats organized, so you do not have to care about hardware specific programming
(e.g. CUDA). We also care about communicating distributed data between processes. On top of this we prepare a set of
basic linear solvers as a jacobi or CG method.

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

 - matrix formats

   - Dense 

   - Sparse

     - CSR

     - COO

     - ELL
     
     - JDS
     
     - DIA

  - backends

    - CPU (OpenMP optimized)

    - GPU

      - CUDA

 - distributions
 
   - blocked

   - cyclic

   - blockcyclic (general block)

   - general
   
 - solvers
 
   - direct solvers
   
     - inverse solver
     
   - iterative solvers
   
     - jacobi method
     
     - conjugated gradiant (CG) method
     
     - successive over-relaxation (SOR) method
     
     - generalized minimal residual (GMRES) method
     
     - simple algebraic multigrid (SAMG) method
    
Work in process that you can find in our feature branches and that will come with next releases are:
 
 - OpenCL-Backend

 - C-Interface

 - PGAS support for communication (till now we support MPI)

 - easy to use configuration of solver through a DSL

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

.. _contact: http://www.libama.org/support.html

Use Cases
=========

So, when is the time to use LAMA?

LAMA is the right decision for you, if you are doing linear algebra on sparse matrices and you want to reach the full
performance of your (parallel) machine without taking care on the kernel code on your own. With LAMA it's easy for you
to write code, that is executable on different heterogeneous machines. 

Possible use cases for LAMA can be found in solving partial elliptic differential equation (e.g. for ground water flow
and oil reservation simulation), image filtering (blur, gauss filtering) and many more. 

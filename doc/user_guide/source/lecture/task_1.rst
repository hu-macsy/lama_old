Task 1: Setup a system of linear equations
==========================================

Linear algebra requires working with objects like scalars, vectors and matrices.
LAMA includes classes, which represent those objects and offers the possibility
to handle with different SparseMatrixFormats ( e.g. CSR, ELL, JDS). This task
should introduce the classes CSRSparseMatrix and DenseVector that have been
already used in task 0 a little bit closer. Additionally the 
SparseAssemblyStorage and some supporting classes are introduced.

Setting up a CSRSparseMatrix
----------------------------

LAMA currently supports the sparse matrix formats CSR, ELLPACK and JDS through
the classes CSRSparseMatrix, ELLSparseMatrix and JDSSparseMatrix. During this
tutorial we will focus on the widely used CSR format.

Lets get started with a good practice: the first essential step is to know how
to create those objects. This Tutorial introduces two possibilities to declare
those objects:

a) To create a CSRSparseMatrix LAMA supports the fileformat \*.mtx from the
Matrix Market [http://math.nist.gov/MatrixMarket/searchtool.html]. This is the
most easiest way to commit data into a Matrix. Download a symetric positive
definite matrix from Matrix Market and create a CSRSparseMatrix from that.

b) If you decide to create a Matrix with your own data you can use a 
SparseAssemblyStorage for this. A SparseAssemblyStorage offers a efficient
possibility to assemble your data in a storage without knowing the quantity of
elements or the order of them. The assembled data can be committed to all
supported sparse matrix formats. To do so you should create an object of type 
SparseAssemblyStorage and fill this storage with elements by using the
set-method. Take care that the SparseMatrix is symetric positive definite so
that the later to implement CG Iterative Solver works. We propose the
discretization of a Possion Equation. To create a CSRSparseMatrix from the 
SparseAssemblyStorage just call the appropriate constructor.


Setting the right hand side and the solution vector
---------------------------------------------------

In the first part of this task has introduced how sparse matrices are handled in
LAMA. For systems of linear equations it is necessary to create vectors as well.
LAMA provides the class DenseVector for this.

Coming back to practice: We have created a filled CSRSparseMatrix. For our
equation system we still need a righthandside-vector and a solution-vector.
Now, your exercise is to create two objects of type DenseVector and fill them
with elements. To fill a DenseVector with data, you have to use a 
HostWriteAccess. The Constructor of this HostWriteAccess expects a LAMA-Array.
You can receive it from your vectors by calling getLocalValues(). With a 
HostWriteAccess it is very simple to set your data using the overloaded
operator[]. The access to data was implemented by the programming technique RAII
(Ressource Acquisition is initialization). This technique is used to assure the
consistency of data across multiple contexts. The concept of a context is
introduced later in this tutorial when it comes to GPU computing. To assure the
consistency we follow the multiple reader single writer idiom, therefore it is 
possible to have multiple ReadAcess, but just one WriteAccess at a time.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_0`", ":doc:`solution_task_1`", ":doc:`task_2`"
:orphan:

Solution Task 1 (a)
===================

.. code-block:: c++
   :emphasize-lines: 9,10

   #include <scai/lama/storage/SparseAssemblyStorage.hpp>
   #include <scai/lama/matrix/CSRSparseMatrix.hpp>
   #include <scai/lama/DenseVector.hpp>

   using namespace lama;

   int main( int argc, char* argv[] ) 
   {
      CSRSparseMatrix<double> m( argv[1] ); /*(1)*/
      IndexType size = m.getNumRows(); /*(2)*/ 
   }

(1) The filename is given as a command-line argument and is the argument of the SparseMatrix-Constructor.
(2) You can get the number of rows by using the method getNumRows().

:download:`Download complete solution Task 1 (a) <../../../solver/examples/lecture/task1a.cpp>`

Solution Task 1 (b)
===================

The following Code is an alternative solution for task 1. Explanations for each
line are listened below.

.. code-block:: c++
   :emphasize-lines: 12,14,29

   #include <scai/lama.hpp>
	
   #include <scai/lama/storage/SparseAssemblyStorage.hpp>
   #include <scai/lama/matrix/CSRSparseMatrix.hpp>
   #include <scai/lama/DenseVector.hpp>

   using namespace lama;

   int main() 
   {
      IndexType size = 4;

      SparseAssemblyStorage<double> sas( size, size, 10 ); /*(1)*/

      for ( IndexType i = 0; i < size; i++ ) /*(2)*/ 
      {
          sas.set( i, i, 2 );
      }
      
      for ( IndexType i = 0; i < size-1; i++ )
      {
         sas.set( i+1, i, 1 );
      }
      
      for ( IndexType i = 0; i < size-1; i++ ) 
      {
         sas.set( i, i+1, 1 );
      }

      CSRSparseMatrix<double> m( sas ); /*(3)*/            

      return 0;
   }

(1) Creation of SparseAssemblyStorage of type double with size 4x4.
(2) Setting some values by using set(). You should only set Non-Zero-Values.
(3) Creation of CSRSparseMatrix of type double and committing the SparseAssemblyStorage.

Setting the right hand side and the solution vector
---------------------------------------------------

.. code-block:: c++
   :emphasize-lines: 10,11,13,18,20

   //#include ...

   using namespace lama;

   int main( int argc, char* argv[] ) 
   {
      CSRSparseMatrix<double> m( argv[1] );
      IndexType size = m.getNumRows();

      DenseVector<double> rhs( size, 0.0 ); /*(1)*/
      HostWriteAccess<double> hwarhs( rhs.getLocalValues() ); /*(2)*/  

      for (int i = 0; i < size; i++ ) /*(3)*/
      {
         hwarhs[i] = i + 1;
      }

      hwarhs.release(); /*(4)*/

      DenseVector<double> solution( size, 0.0 ); /*(5)*/
    }

(1) Creation of DenseVector rhs of type double and default-values 0.0.
(2) Creation of HostWriteAccess of type double for DenseVector rhs. The Constructor requires a LAMA-Array. You can get it by calling the getLocalValues()-method of your DenseVector.
(3) Setting values of rhs by yourself. The overloaded operator[] makes it easy to handle it.
(4) Release of HostWriteAccesses. Instead of releasing the HostWriteAccess you can use a block { /\* set() here \*/ }. The release()-method will be automatically called of the Destructor at the end of this block.
(5) Creation of DenseVector solution. Default-value is 0.0.

:download:`Download complete solution Task 1 (b) <../../../solver/examples/lecture/task1b.cpp>`

.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_1`", ":doc:`../lecture`", ":doc:`task_2`"

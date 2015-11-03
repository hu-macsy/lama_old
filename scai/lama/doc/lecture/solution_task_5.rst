:orphan:

Solution Task 5
===============

Setting a Context was realized as simple as possible:

.. code-block:: c++
   :emphasize-lines: 8,16

   //#include {...}

   using namespace lama;

   int main( int argc, char* argv[] )
   {
   
      lama::ContextPtr cudaContext = ContextFactory::getContext( Context::CUDA, 0 ); /*(1)*/     

      CSRSparseMatrix<double> m( argv[1] );
      IndexType size = m.getNumRows();
   
      //Declaration of objects
      {...} 

      m.setContext( cudaContext ); /*(2)*/
      rhs.setContext( cudaContext );
      solution.setContext( cudaContext );

      {...}

      return 0;
   }

(1) : Getting a CudaContext for cuda device 0.
(2) : Setting a CudaContext to matrix and vectors.

:download:`Download complete solution Task 5 <../../../lama/examples/lecture/task5.cpp>`

.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_5`", ":doc:`../lecture`", ":doc:`task_6`"

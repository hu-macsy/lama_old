.. _tutorial_solution_task5:

Solution of Task 5
==================

Setting a Context was realized as simple as possible:

.. code-block:: c++
   :emphasize-lines: 8,16

   //#include {...}

   using namespace lama;

   int main( int argc, char* argv[] )
   {
   
      /*(1)*/     lama::ContextPtr cudaContext = ContextFactory::getContext( Context::CUDA, 0 ); 

      CSRSparseMatrix<double> m( argv[1] );
      IndexType size = m.getNumRows();
   
      //Declaration of objects
      {...} 

      /*(2)*/     m.setContext( cudaContext );
      rhs.setContext( cudaContext );
      solution.setContext( cudaContext );

      {...}

      return 0;
   }

(1) : Getting a CudaContext for cuda device 0.
(2) : Setting a CudaContext to matrix and vectors.

The solution of task 5 can be downloaded `here`__ 

__ http://libama.sourceforge.net/tutorial/solutions/task5.cpp


.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_5`", ":doc:`index`", ":doc:`task_6`"

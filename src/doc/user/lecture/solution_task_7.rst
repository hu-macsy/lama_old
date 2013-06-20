:orphan:

Solution Task 7
===============

A possible solution for task 7 is shown in the following sourcecode:

::

   #include <...>

   using namespace lama;

   int main()
   {

   LAMA_REGION( "CUDA-Region" );
   lama::ContextPtr mCudaContext = ContextFactory::getContext( Context::CUDA, 0 );

   {...}

   //Creating CG-Solver

   LAMA_REGION_START( "CG-Solver-Region" );
   for (int k = 1; k < 10; k++) 
   {
      if (norm(r) < eps) 
      {
         break;
      }
      else
      {
         Ad = m * d;
         alpha = rOld / ( d * Ad );

         solution = solution + ( alpha * d );
         r = r - alpha * Ad;

         rNew = r * r;
         beta = rNew / rOld;

         d = r + beta * d;
         rOld = rNew;
      }
   }
   LAMA_REGION_END( "CG-Solver-Region" );

   {...}
   }

The performance can be analyzed with Vampir. Especially if you run your program
MPI parallel, the communication between your processes will be presented.

Starting Vampir:

.. code-block:: bash

   vampir foo.otf
   
.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_7`", ":doc:`../lecture`", "-"
   
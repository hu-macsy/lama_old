.. _tutorial_solution_task7:

Solution of Task 7
==================

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

Starting Vampir on NOVA cluster:

.. code-block:: bash

   module load UNITE
   module load vampir

Open a new bash on your local machine and connect to NOVA via:

.. code-block:: bash

   ssh -X <user>@129.184.111.23
   ssh -X nv47
   echo $DISPLAY

Switch back to the original qsub-shell and export the value - e.g.:

.. code-block:: bash
 
   export DISPLAY=localhost:13.0

Start vampir in your qsub-shell:

.. code-block:: bash
   
   vampir &

   
.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_7`", ":doc:`index`", "-"
   
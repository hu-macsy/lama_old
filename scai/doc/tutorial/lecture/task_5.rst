:orphan:

Task 5: Let the CG Iterative Solver run on a GPU
------------------------------------------------

LAMA provides the concept of contexts that allows to run code on different compute devices.
A LAMA-Context represents a compute and
data device, which can be set by a developer to let the code run on this device.
For GPU devices this context is called CUDA context as it uses for the implementation of
kernel routines the CUDA language.

Back to business: This task should exemplarily demonstrate how to run a
CG-Solver on a GPU with CUDA. For this we need to create a CUDA context.
To do so we can obtain a ContextPtr from the ContextFactory, by querying for a
CUDA context. 

.. code-block:: c++

    ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA, 0 );

For more detailed explanation of this, refer to :ref:`scaihmemo:hmemo-contextFactory`.

If we set this context to be the context for all matrices and
vectors involved in the CG-Solver it implicitly runs on the GPU. You can set the
context of a matrix or a vector either by passing it in a constructor or
by calling their ``setContextPtr`` method.

**Excercise**: Run the CG solver code on a GPU.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_4`", ":doc:`solution_task_5`", ":doc:`task_6`"

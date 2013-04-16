:orphan:

Task 5: Let the CG Iterative Solver run on a GPU
================================================

LAMA provides the concept of contextes to run code on different compute devices.
For GPU, LAMA supports OpenCL and CUDA. A LAMA-Context represents a compute and
data device, which can be set by a developer to let the code run on this device.

Back to business: This task should exemplarily demonstrate how to run your
developed CG-Solver on GPU with CUDA. For this we need to create a CUDA context.
To do so we can obtain an ContextPtr from the ContextFactory, by querying for a
CUDA context. If we set this context to be the context for all matrices and
vectors involved in the CG-Solver it implicitly runs on the GPU. You can set the
context of a matrix or a vector by calling their setContext method.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_4`", ":doc:`solution_task_5`", ":doc:`task_6`"
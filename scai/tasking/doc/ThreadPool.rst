ThreadPool
==========

 * Class to define a pool of threads
 * Avoids overhead of thread creation for asynchronous execution
 * One thread pool is used for TaskSyncToken 

.. code-block:: c++

   ThreadPool pool( 10 );


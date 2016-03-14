ThreadPool
==========

 * Class to define a pool of threads
 * Avoids overhead of thread creation for asynchronous execution

.. code-block:: c++

    ThreadPool pool( n );   // creates a pool with n threads

    std::cout << "ThreadPool, size = " << pool.size() << " created" << std::endl;

    // schedule 3 tasks, task3 must wait for completion of task1 or task2 

    shared_ptr<ThreadPoolTask> task1 = pool.schedule( bind( &work, 1 ) );
    shared_ptr<ThreadPoolTask> task2 = pool.schedule( bind( &work, 2 ) );
    shared_ptr<ThreadPoolTask> task3 = pool.schedule( bind( &work, 3 ) );

    // wait for completion

    pool.wait( task1 );
    pool.wait( task2 );
    pool.wait( task3 );

Some remarks:

 * Functions with binded arguments can be scheduled within a thread pool 
 * A struct of type ``ThreadPoolTask`` is returned that contains information about
   the state of the scheduled function (e.g. SCHEDULED, RUNNING, FINISHED)
 * The struct can be used to wait for its completion
 * The struct of ``ThreadPoolTask`` is always returned via a shared pointer as 
   ownership can be taken either by the pool (until its finished) or by the application.

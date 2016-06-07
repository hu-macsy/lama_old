.. _Task:

Task
====

 * Encapsulates a function call scheduled in a thread pool
 * Constructor schedules the function
 * Destructor wais for its completion.

.. code-block:: c++

    using namespace scai;

    {
        tasking::Task t1( common::bind( &work, 3 ) );
        tasking::Task t2( common::bind( &work, 2 ) );
        tasking::Task t3( common::bind( &work, 1 ) );
        ...
        // implicit wait at the end of the scope
    }
    // here all tasks are finished

Some remarks:

 * Size of the thread pool can be configured by the environment variable ``SCAI_THREADPOOL_SIZE``

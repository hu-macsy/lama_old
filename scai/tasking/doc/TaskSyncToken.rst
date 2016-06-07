.. _TaskSyncToken:

TaskSyncToken
=============

 * Derived class of SyncToken 
 * Stands for the execution of a function by a separate thread.
 * Functionality is very close to the class Task

Here is a typical example of its use:

.. code-block:: c++

   void function f( ... )
   {
       ....
   }

   ...
   {
        common::unique_ptr<SyncToken> token( new TaskSyncToken( common::bind( &f, arg1, .... ) ) );
        ....
        // implicit synchronization at end of the scope
   }



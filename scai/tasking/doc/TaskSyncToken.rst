TaskSyncToken
=============

 * Derived class of SyncToken 
 * Stands for the execution of a function by a separate thread.

Here is a typical example of its use:

.. code-block:: c++

   void function f()
   {
       ....
   }

   ...
   {
        common::unique_ptr<SyncToken> token( new TaskSyncToken( f ) );
        ....
        // implicit synchronization at end of the scope
   }



SyncToken
=========

 * Base class
 * One derived class for each kind of asynchronous execution

Here is a typical example of its use:

.. code-block:: c++

   SyncToken* token = new DerivedSyncToken( .... );  // starts asynchronous execution
   ...   // here other work can be done
   token->wait();  // sychronize, i.e. wait for completion of the asychronous task


The wait operation can be called explicitly, but the synchronization is also done
implicitly with the destructor of the token.

It is recommended to use a smart pointer for the token to guarantee that the task is 
really completed.

.. code-block:: c++

   {
        common::unique_ptr<SyncToken> token( new DerivedSyncToken( .... ) );
        ....
        // implicit synchronization at end of the scope
   }

Very important: it is possible to add finalize routines for the token:

.. code-block:: c++

    token->pushRoutine( ... );   // will be called at wait after the async op is finished.



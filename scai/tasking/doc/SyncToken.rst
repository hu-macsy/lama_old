SyncToken
=========

 * Base class
 * One derived class for each kind of asynchronous execution

Here is a typical example of its use:

.. code-block:: c++

   // constructor of a derived SyncToken class starts asynchronous execution

   SyncToken* token = new DerivedSyncToken( .... );  

   ...   // here other work can be done

   token->wait();  // sychronize, i.e. wait for completion of the asychronous task

The wait operation can be called explicitly, but the synchronization is also done
implicitly by the destructor of the token.

.. code-block:: c++

   SyncToken* token = new DerivedSyncToken( .... );  

   ...

   delete token;   // implicit synchronization

It is recommended to use a smart pointer for the token to guarantee that the
destructor of the object is called and the synchronization is really done.
Furthermore, the SyncToken object might call other cleanup routines whose
calls are mandatory.

.. code-block:: c++

   {
        common::unique_ptr<SyncToken> token( new DerivedSyncToken( .... ) );

        ....

        // implicit synchronization at end of the scope with destructor
   }

At the end of an asynchronous operation it is often necessary to release allocated
resources, e.g. memory or accesses to devices. Threrefore it is possible to add
function calls to a SyncToken where these function calls will be executed after the synchronization.
This operation is usually called by the constructor of derived classes.

.. code-block:: c++

    void freeResources( T1 arg1, .., Tn argn )
    {
       ...
    }

    ...

    common::function<void()> f = common::bind( freeResources, val1, .., valn );
    token->pushRoutine( f );   // will be called at wait after the async op is finished.


.. _NoSyncToken:

NoSyncToken
===========

 * Derived class for synchronous execution

There might be virtual routines in classes that return a SyncToken for
an asynchronous execution.

.. code-block:: c++

    SyncToken* memcopy( ... )

Derived classes might redefine such a routine and return a SyncToken that 
fits best. The class NoSyncToken is an alternative for all situations

* where the overhead is too high to run an asynchronous execution
* where an asynchronous execution is not supported or not possible

Then the above function might return a NoSyncToken that is just a dummy
object that provides also the wait operation. Actually the execution itself
is already finished before the above function returns.
